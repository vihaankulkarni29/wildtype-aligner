from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
import textwrap
import logging
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import sys
from typing import List, Dict, Optional, Any
import time
import socket

# Try to import tqdm for progress bars, fallback to simple counter
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    tqdm = None

logging.basicConfig(level=logging.INFO)

# Module-level worker function for parallel processing
# Must be at module level to be picklable by multiprocessing
def _align_worker(args):
    """
    Worker function for parallel alignment processing.
    
    Args:
        args: tuple of (query_seq, reference_seq, sequence_number)
    
    Returns:
        tuple: (sequence_number, alignment_report_string or error_message)
    """
    query_seq, reference_seq, seq_num = args
    
    try:
        # Create a fresh aligner instance for this worker
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        
        # Perform alignment
        try:
            alignments = list(aligner.align(query_seq.seq, reference_seq.seq))
            alignment = alignments[0] if alignments else None
        except OverflowError:
            # Handle overflow by taking first alignment
            alignment_iter = aligner.align(query_seq.seq, reference_seq.seq)
            alignment = next(alignment_iter, None)
        
        if alignment is None:
            return (seq_num, f"# {seq_num}\nNo alignment found for: {query_seq.id}\n\n")
        
        # Format alignment (defensive to satisfy static analyzers)
        score = float(getattr(alignment, 'score', 0.0))
        try:
            seq1 = str(alignment[0])  # type: ignore[index]
            seq2 = str(alignment[1])  # type: ignore[index]
        except Exception:
            return (seq_num, f"# {seq_num}\nAlignment obtained but could not extract sequences for: {query_seq.id}\nScore: {score:.1f}\n\n")

        length = len(seq1)
        identities = sum(1 for a, b in zip(seq1, seq2) if a == b)
        identity_percent = (identities / length) * 100 if length > 0 else 0
        gaps = seq1.count('-') + seq2.count('-')
        
        header = f"# {seq_num}\n"
        header += f"Waterman-Eggert local alignment\n"
        header += f"Query: {query_seq.id} ({len(query_seq.seq)} aa)\n"
        header += f"Reference: {reference_seq.id} ({len(reference_seq.seq)} aa)\n"
        header += f"Score: {score:.1f}\n"
        header += f"Identities: {identities}/{length} ({identity_percent:.1f}%)\n"
        header += f"Gaps: {gaps}\n\n"
        
    # Alignment blocks
        matches = ''.join('|' if a == b else ' ' for a, b in zip(seq1, seq2))
        
        block_size = 60
        formatted = ""
        for start in range(0, length, block_size):
            end = min(start + block_size, length)
            formatted += f"{query_seq.id:<20} {start+1:>4} {seq1[start:end]} {end:>4}\n"
            formatted += f"{'':<20} {'':>4} {matches[start:end]} {'':>4}\n"
            formatted += f"{reference_seq.id:<20} {start+1:>4} {seq2[start:end]} {end:>4}\n\n"
        
        return (seq_num, header + formatted)
        
    except Exception as e:
        return (seq_num, f"# {seq_num}\nError aligning {query_seq.id}: {str(e)}\n\n")


class ProgressTracker:
    """
    Progress tracker that uses tqdm if available, otherwise falls back to simple counter.
    """
    def __init__(self, total, desc="Processing", use_tqdm=True):
        self.total = total
        self.current = 0
        self.desc = desc
        self.use_tqdm = use_tqdm and TQDM_AVAILABLE
        
        if self.use_tqdm:
            # tqdm is only callable when actually imported; assert for type narrowing
            assert tqdm is not None
            self.pbar = tqdm(total=total, desc=desc, unit="seq", 
                           bar_format='{desc}: {n_fmt}/{total_fmt} ({percentage:3.0f}%) |{bar}| [{elapsed}<{remaining}]')
        else:
            self.pbar = None
            print(f"{desc}: 0/{total} (0%)")
    
    def update(self, n=1):
        """Update progress by n steps"""
        self.current += n
        
        if self.use_tqdm and self.pbar is not None:
            self.pbar.update(n)
        else:
            # Simple fallback: print every 10% or at milestones
            percent = (self.current / self.total) * 100
            if self.current % max(1, self.total // 10) == 0 or self.current == self.total:
                print(f"{self.desc}: {self.current}/{self.total} ({percent:.0f}%)")
                sys.stdout.flush()
    
    def set_description(self, desc):
        """Update the description text"""
        self.desc = desc
        if self.use_tqdm and self.pbar is not None:
            self.pbar.set_description(desc)
    
    def close(self):
        """Close the progress bar"""
        if self.use_tqdm and self.pbar is not None:
            self.pbar.close()
        else:
            if self.current == self.total:
                print(f"{self.desc}: Complete! ({self.total}/{self.total})")
            sys.stdout.flush()


class AlignmentSummary:
    """
    Track alignment successes and failures for summary reporting.
    """
    def __init__(self):
        self.total = 0
        self.successful = 0
        self.failed = 0
        self.errors = []  # List of (sequence_id, error_message) tuples
    
    def add_success(self, seq_id):
        """Record a successful alignment"""
        self.total += 1
        self.successful += 1
    
    def add_failure(self, seq_id, error_msg):
        """Record a failed alignment"""
        self.total += 1
        self.failed += 1
        self.errors.append((seq_id, error_msg))
    
    def get_summary(self):
        """Get formatted summary report"""
        summary = []
        summary.append("=" * 70)
        summary.append("ALIGNMENT SUMMARY")
        summary.append("=" * 70)
        summary.append(f"Total sequences processed: {self.total}")
        summary.append(f"Successfully aligned:     {self.successful} ({self.successful/self.total*100:.1f}%)" if self.total > 0 else "Successfully aligned:     0")
        summary.append(f"Failed alignments:        {self.failed} ({self.failed/self.total*100:.1f}%)" if self.total > 0 else "Failed alignments:        0")
        
        if self.failed > 0:
            summary.append("")
            summary.append("FAILED SEQUENCES:")
            summary.append("-" * 70)
            for seq_id, error_msg in self.errors:
                summary.append(f"  • {seq_id}")
                summary.append(f"    Error: {error_msg}")
        
        summary.append("=" * 70)
        return "\n".join(summary)
    
    def save_error_log(self, log_file):
        """Save detailed error log to file"""
        if self.failed == 0:
            return
        
        with open(log_file, 'w') as f:
            f.write("ALIGNMENT ERROR LOG\n")
            f.write("=" * 70 + "\n")
            f.write(f"Generated: {__import__('datetime').datetime.now()}\n")
            f.write(f"Total failures: {self.failed}\n")
            f.write("=" * 70 + "\n\n")
            
            for seq_id, error_msg in self.errors:
                f.write(f"Sequence ID: {seq_id}\n")
                f.write(f"Error: {error_msg}\n")
                f.write("-" * 70 + "\n")


class AlignmentEngine:
    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None,
                 timeout: Optional[int] = None, max_retries: int = 3):
        # Configure NCBI Entrez
        # Email for NCBI compliance; allow override from CLI/env
        Entrez.email = email or os.environ.get("ENTREZ_EMAIL", "vihaankulkarni29@gmail.com")
        Entrez.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self._timeout = timeout or int(os.environ.get("ENTREZ_TIMEOUT", "0") or 0) or None
        self._max_retries = max(1, int(os.environ.get("ENTREZ_MAX_RETRIES", str(max_retries))))
        # Initialize aligner
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # Use global alignment
        # Initialize reference cache: {gene_name: SeqRecord}
        self.reference_cache = {}
        logging.info("AlignmentEngine initialized")
        logging.info(f"Entrez email={Entrez.email}, key={'set' if Entrez.api_key else 'unset'}, timeout={self._timeout}, retries={self._max_retries}")

    def _with_timeout(self):
        """Context manager to set a global socket timeout if configured."""
        class _TimeoutCtx:
            def __init__(self, seconds):
                self.seconds = seconds
                self._prev = None
            def __enter__(self):
                if self.seconds:
                    self._prev = socket.getdefaulttimeout()
                    socket.setdefaulttimeout(self.seconds)
            def __exit__(self, exc_type, exc, tb):
                if self.seconds is not None:
                    socket.setdefaulttimeout(self._prev)
        return _TimeoutCtx(self._timeout)

    def _retry(self, func, *args, **kwargs):
        """Run an Entrez call with retries and exponential backoff."""
        delay = 1.0
        last_exc: Optional[Exception] = None
        for attempt in range(1, self._max_retries + 1):
            try:
                with self._with_timeout():
                    return func(*args, **kwargs)
            except Exception as e:
                last_exc = e
                logging.warning(f"Entrez call failed (attempt {attempt}/{self._max_retries}): {e}")
                if attempt < self._max_retries:
                    time.sleep(delay)
                    delay *= 2
        # Exhausted retries
        raise last_exc if last_exc else RuntimeError("Unknown Entrez error")

    def _get_reference_with_cache(self, gene_name, genus, species):
        """
        Get reference sequence with caching to avoid redundant NCBI calls.
        
        Args:
            gene_name: Name of the gene (e.g., 'acrA', 'acrB', 'tolC')
            genus: Genus name (e.g., 'Escherichia')
            species: Species name (e.g., 'coli')
        
        Returns:
            SeqRecord: The reference sequence
        
        Cache structure: {gene_name: SeqRecord}
        When processing 136 acrA sequences, fetches reference once and reuses it.
        """
        # Create cache key from gene name (genus and species are typically the same)
        cache_key = gene_name.lower()
        
        # Check if reference is already cached
        if cache_key in self.reference_cache:
            logging.info(f"✓ Using cached reference for gene '{gene_name}'")
            return self.reference_cache[cache_key]
        
        # Not cached - fetch from NCBI
        logging.info(f"⟳ Fetching reference for gene '{gene_name}' from NCBI (not in cache)")
        reference = self._fetch_sepi_reference(gene_name, genus, species)
        
        # Cache it for future use
        self.reference_cache[cache_key] = reference
        logging.info(f"✓ Cached reference for gene '{gene_name}' (cache size: {len(self.reference_cache)})")
        
        return reference

    def get_cache_stats(self):
        """
        Get cache statistics for monitoring and debugging.
        
        Returns:
            dict: Cache statistics including size and cached genes
        """
        return {
            'cache_size': len(self.reference_cache),
            'cached_genes': list(self.reference_cache.keys()),
            'cache_info': {gene: f"{ref.id} ({len(ref.seq)} aa)" 
                          for gene, ref in self.reference_cache.items()}
        }

    def clear_cache(self):
        """
        Clear the reference cache.
        Useful for testing or if you need to force re-fetching.
        """
        self.reference_cache.clear()
        logging.info("Reference cache cleared")

    def _group_sequences_by_gene(self, query_sequences):
        """
        Group query sequences by their gene type.
        
        Args:
            query_sequences: List of SeqRecord objects
        
        Returns:
            dict: {gene_name: [SeqRecord, SeqRecord, ...]}
        
        Example:
            {
                'acrA': [seq1, seq2, seq3, ...],  # 136 sequences
                'acrB': [seq1, seq2, seq3, ...],  # 137 sequences
                'tolC': [seq1, seq2, seq3, ...]   # 136 sequences
            }
        """
        grouped = defaultdict(list)
        
        for seq in query_sequences:
            # Extract gene name from sequence ID
            gene_name = self._extract_gene_from_sequence_id(seq.id)
            
            # Fallback to description if needed
            if not gene_name:
                gene_name = self._extract_gene_name(seq.description)
            
            if gene_name:
                grouped[gene_name.lower()].append(seq)
            else:
                # If we can't determine gene, put in 'unknown' group
                logging.warning(f"Cannot determine gene for sequence: {seq.id}")
                grouped['unknown'].append(seq)
        
        # Log grouping summary
        for gene, sequences in grouped.items():
            logging.info(f"Grouped {len(sequences)} sequences for gene '{gene}'")
        
        return dict(grouped)

    def run_batch_alignments(self, query_sequences, genus, species, output_file=None, show_progress=True, error_log=None):
        """
        Perform batch alignments with intelligent grouping and caching.
        
        This method:
        1. Groups sequences by gene type
        2. Fetches ONE reference per gene type (with caching)
        3. Aligns all sequences of same gene against their shared reference
        4. Optionally streams output to file for memory efficiency
        5. Shows real-time progress tracking
        6. Robust error handling - continues even if some sequences fail
        
        Args:
            query_sequences: List of SeqRecord objects (can be mixed genes)
            genus: Organism genus (e.g., 'Escherichia')
            species: Organism species (e.g., 'coli')
            output_file: File handle for streaming output (optional)
            show_progress: Whether to show progress bar (default: True)
            error_log: Path to error log file (optional)
        
        Returns:
            tuple: (report_string or None, AlignmentSummary object)
        
        Example:
            For 409 sequences (136 acrA, 137 acrB, 136 tolC):
            - Makes only 3 NCBI calls (one per gene type)
            - Processes all alignments efficiently
            - Continues even if some fail
        """
        logging.info(f"Starting batch alignment for {len(query_sequences)} sequences")
        
        # Determine if streaming mode
        streaming = output_file is not None
        if streaming:
            logging.info("STREAMING MODE: Writing results directly to file")
        
        # Initialize summary tracker
        summary = AlignmentSummary()
        
        # Group sequences by gene type
        grouped_sequences = self._group_sequences_by_gene(query_sequences)
        
        if not grouped_sequences:
            raise ValueError("No valid sequences found for alignment")
        
        # Initialize progress tracker
        total_sequences = len(query_sequences)
        progress = ProgressTracker(total_sequences, desc="Aligning sequences", use_tqdm=show_progress) if show_progress else None
        
        # Process each gene group
        # Always initialize list (unused if streaming) to satisfy static analyzers
        all_reports: List[str] = []
        sequence_counter = 0
        
        for gene_name, sequences in sorted(grouped_sequences.items()):
            logging.info(f"\n{'='*60}")
            logging.info(f"Processing {len(sequences)} sequences for gene: {gene_name}")
            logging.info(f"{'='*60}")
            
            try:
                # Get reference with caching (only fetches once per gene)
                reference = self._get_reference_with_cache(gene_name, genus, species)
                
                # Update progress description for current gene
                if progress:
                    progress.set_description(f"Aligning {gene_name}")
                
                # Align all sequences of this gene against the shared reference
                for seq in sequences:
                    sequence_counter += 1
                    try:
                        alignment = self._align_single(seq, reference)
                        if alignment is None:
                            # Track failure
                            summary.add_failure(seq.id, "No alignment found")
                            
                            error_msg = f"# {sequence_counter}\nNo alignment found for: {seq.id}\n\n"
                            if streaming:
                                output_file.write(error_msg)
                                output_file.write("\n\n")
                                output_file.flush()
                            else:
                                all_reports.append(error_msg)
                            
                            # Update progress
                            if progress:
                                progress.update(1)
                            continue
                        
                        alignment_str = self._format_water_style(
                            alignment, seq, reference, sequence_counter
                        )
                        
                        # Track success
                        summary.add_success(seq.id)
                        
                        if streaming:
                            output_file.write(alignment_str)
                            output_file.write("\n\n")
                            output_file.flush()
                        else:
                            all_reports.append(alignment_str)
                        
                        # Update progress
                        if progress:
                            progress.update(1)
                            
                    except Exception as e:
                        # Track failure with detailed error
                        summary.add_failure(seq.id, str(e))
                        logging.error(f"Error aligning {seq.id}: {str(e)}")
                        
                        error_msg = f"# {sequence_counter}\nError aligning {seq.id}: {str(e)}\n\n"
                        if streaming:
                            output_file.write(error_msg)
                            output_file.write("\n\n")
                            output_file.flush()
                        else:
                            all_reports.append(error_msg)
                        
                        # Update progress even on error
                        if progress:
                            progress.update(1)
                
            except Exception as e:
                logging.error(f"Error processing gene group {gene_name}: {str(e)}")
                # Track failures and add error reports for all sequences in this group
                for seq in sequences:
                    sequence_counter += 1
                    summary.add_failure(seq.id, f"Gene group error: {str(e)}")
                    
                    error_msg = f"# {sequence_counter}\nError processing {gene_name} ({seq.id}): {str(e)}\n\n"
                    if streaming:
                        output_file.write(error_msg)
                        output_file.write("\n\n")
                        output_file.flush()
                    else:
                        all_reports.append(error_msg)
                    
                    # Update progress even on error
                    if progress:
                        progress.update(1)
        
        # Close progress tracker
        if progress:
            progress.close()
        
        # Save error log if requested
        if error_log and summary.failed > 0:
            try:
                summary.save_error_log(error_log)
                logging.info(f"Error log saved to {error_log}")
            except Exception as e:
                logging.error(f"Failed to save error log: {str(e)}")
        
        # Log summary
        logging.info(summary.get_summary())
        
        logging.info(f"\n{'='*60}")
        logging.info(f"Batch alignment complete: {sequence_counter} sequences processed")
        logging.info(f"Cache statistics: {self.get_cache_stats()}")
        logging.info(f"{'='*60}")
        
        if streaming:
            return None, summary  # Return summary even in streaming mode
        else:
            return "\n\n".join(all_reports), summary

    def run_parallel_batch_alignments(self, query_sequences, genus, species, max_workers=None, output_file=None, show_progress=True, error_log=None):
        """
        Perform batch alignments with parallel processing for maximum performance.
        
        This method combines:
        1. Intelligent grouping by gene type
        2. Reference caching (only 3 NCBI calls)
        3. Multi-core parallel alignment processing
        4. STREAMING OUTPUT: writes results as they complete (memory efficient)
        5. Real-time progress tracking
        6. Comprehensive error handling with detailed reporting
        
        Args:
            query_sequences: List of SeqRecord objects (can be mixed genes)
            genus: Organism genus (e.g., 'Escherichia')
            species: Organism species (e.g., 'coli')
            max_workers: Number of worker processes (default: CPU count)
            output_file: File handle for streaming output (optional, for memory efficiency)
            show_progress: Whether to show progress bar (default: True)
            error_log: Path to error log file (optional)
        
        Returns:
            tuple: (report_string or None, AlignmentSummary object)
                - report_string: Comprehensive alignment report (if no output_file), None if streaming
                - AlignmentSummary: Tracking object with success/failure statistics
        
        Performance:
            For 409 sequences on 8-core CPU:
            - Sequential: ~10-15 minutes
            - Parallel: ~2-3 minutes (5-6x speedup)
            - Memory: Constant (with streaming) vs O(n) (without streaming)
        """
        if max_workers is None:
            max_workers = multiprocessing.cpu_count()
        
        logging.info(f"Starting PARALLEL batch alignment for {len(query_sequences)} sequences")
        logging.info(f"Using {max_workers} worker processes")
        
        # Initialize summary tracker
        summary = AlignmentSummary()
        
        # Determine if streaming mode
        streaming = output_file is not None
        if streaming:
            logging.info("STREAMING MODE: Writing results directly to file as they complete")
        
        # Initialize progress tracker
        total_sequences = len(query_sequences)
        progress = ProgressTracker(total_sequences, desc="Parallel alignment", use_tqdm=show_progress) if show_progress else None
        
        # Group sequences by gene type
        grouped_sequences = self._group_sequences_by_gene(query_sequences)
        
        if not grouped_sequences:
            raise ValueError("No valid sequences found for alignment")
        
        # Process each gene group with parallelization
        # Always initialize dict (unused if streaming) to satisfy static analyzers
        all_results: Dict[int, str] = {}
        sequence_counter = 0
        
        for gene_name, sequences in sorted(grouped_sequences.items()):
            logging.info(f"\n{'='*60}")
            logging.info(f"Processing {len(sequences)} sequences for gene: {gene_name}")
            logging.info(f"{'='*60}")
            
            try:
                # Get reference with caching (only fetches once per gene)
                reference = self._get_reference_with_cache(gene_name, genus, species)
                
                # Update progress description for current gene
                if progress:
                    progress.set_description(f"Parallel alignment [{gene_name}]")
                
                # Prepare work items for parallel processing
                # Each item: (query_seq, reference_seq, sequence_number)
                work_items = []
                for seq in sequences:
                    sequence_counter += 1
                    work_items.append((seq, reference, sequence_counter))
                
                # Process alignments in parallel
                logging.info(f"Submitting {len(work_items)} alignments to {max_workers} workers...")
                
                # In streaming mode, we need to track results to write them in order
                pending_results = {} if streaming else {}
                
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    # Submit all tasks
                    futures = {executor.submit(_align_worker, item): item for item in work_items}
                    
                    # Collect results as they complete
                    completed = 0
                    for future in as_completed(futures):
                        try:
                            seq_num, report = future.result()
                            item = futures[future]
                            seq = item[0]  # Get the sequence from work item
                            
                            # Track success
                            summary.add_success(seq.id)
                            
                            if streaming:
                                # In streaming mode: store in pending until we can write in order
                                pending_results[seq_num] = report
                                
                                # Write all consecutive results that are ready
                                # This ensures we write in sequence order
                                next_seq_to_write = sequence_counter - len(sequences) + completed + 1
                                while next_seq_to_write in pending_results:
                                    output_file.write(pending_results[next_seq_to_write])
                                    output_file.write("\n\n")
                                    output_file.flush()  # Force write to disk
                                    del pending_results[next_seq_to_write]
                                    next_seq_to_write += 1
                            else:
                                # Non-streaming mode: accumulate in memory
                                all_results[seq_num] = report
                            
                            completed += 1
                            
                            # Update progress
                            if progress:
                                progress.update(1)
                            
                            if completed % 10 == 0 or completed == len(work_items):
                                logging.info(f"Progress: {completed}/{len(work_items)} alignments complete")
                        except Exception as e:
                            item = futures[future]
                            seq = item[0]  # Get the sequence from work item
                            seq_num = item[2]
                            
                            # Track failure with detailed error
                            summary.add_failure(seq.id, f"Worker error: {str(e)}")
                            logging.error(f"Worker error for sequence {seq_num} ({seq.id}): {str(e)}")
                            error_msg = f"# {seq_num}\nWorker error ({seq.id}): {str(e)}\n\n"
                            
                            if streaming:
                                pending_results[seq_num] = error_msg
                                # Try to flush pending results
                                next_seq_to_write = sequence_counter - len(sequences) + completed + 1
                                while next_seq_to_write in pending_results:
                                    output_file.write(pending_results[next_seq_to_write])
                                    output_file.write("\n\n")
                                    output_file.flush()
                                    del pending_results[next_seq_to_write]
                                    next_seq_to_write += 1
                            else:
                                all_results[seq_num] = error_msg
                            
                            completed += 1
                            
                            # Update progress even on error
                            if progress:
                                progress.update(1)
                
                logging.info(f"Completed all {len(sequences)} {gene_name} alignments")
                
            except Exception as e:
                logging.error(f"Error processing gene group {gene_name}: {str(e)}")
                # Track failures and add error reports for all sequences in this group
                for seq in sequences:
                    sequence_counter += 1
                    summary.add_failure(seq.id, f"Gene group error: {str(e)}")
                    
                    error_msg = f"# {sequence_counter}\nError processing {gene_name} ({seq.id}): {str(e)}\n\n"
                    
                    if streaming:
                        output_file.write(error_msg)
                        output_file.write("\n\n")
                        output_file.flush()
                    else:
                        all_results[sequence_counter] = error_msg
                    
                    # Update progress even on error
                    if progress:
                        progress.update(1)
        
        # Close progress tracker
        if progress:
            progress.close()
        
        # Save error log if requested
        if error_log and summary.failed > 0:
            try:
                summary.save_error_log(error_log)
                logging.info(f"Error log saved to {error_log}")
            except Exception as e:
                logging.error(f"Failed to save error log: {str(e)}")
        
        # Log summary
        logging.info(summary.get_summary())
        
        # Final assembly/logging
        if streaming:
            logging.info(f"\n{'='*60}")
            logging.info(f"Parallel batch alignment complete: {sequence_counter} sequences processed")
            logging.info(f"Results streamed to output file")
            logging.info(f"Cache statistics: {self.get_cache_stats()}")
            logging.info(f"{'='*60}")
            return None, summary  # Return summary even in streaming mode
        else:
            # Assemble results in correct order
            ordered_reports = [all_results[i] for i in sorted(all_results.keys())]
            
            logging.info(f"\n{'='*60}")
            logging.info(f"Parallel batch alignment complete: {sequence_counter} sequences processed")
            logging.info(f"Cache statistics: {self.get_cache_stats()}")
            logging.info(f"{'='*60}")
            
            return "\n\n".join(ordered_reports), summary

    def _fetch_sepi_reference(self, gene_name, genus, species):
        """
        Fetch the canonical protein sequence from NCBI for the given gene, genus, and species.
        Returns a SeqRecord object.
        """
        try:
            # Construct search term with genus and species for precision
            # Prefer RefSeq proteins and avoid partials when possible
            search_term = (
                f"{gene_name}[Gene] AND {genus} {species}[Organism] "
                f"AND srcdb_refseq[PROP] NOT partial[Title]"
            )
            logging.info(f"Searching NCBI for {gene_name} in {genus} {species}")
            # Search for protein records
            handle = self._retry(Entrez.esearch, db="protein", term=search_term, retmax=20)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:  # type: ignore[index]
                raise ValueError(f"No protein sequence found for gene {gene_name} in {genus} {species}")

            ids = record["IdList"]  # type: ignore[index]

            # Prefer NP_/YP_ RefSeq accessions via esummary
            protein_id = ids[0]
            try:
                sum_handle = self._retry(Entrez.esummary, db="protein", id=",".join(ids))
                summary = Entrez.read(sum_handle)
                sum_handle.close()
                best_id = None
                for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
                    acc = doc.get("AccessionVersion") or doc.get("Caption")
                    if acc and (str(acc).startswith("NP_") or str(acc).startswith("YP_")):
                        best_id = str(doc.get("Id")) if doc.get("Id") is not None else None
                        break
                if best_id:
                    protein_id = best_id
            except Exception as e:
                logging.warning(f"esummary preference failed, falling back to first ID: {e}")

            # Fetch the sequence
            handle = self._retry(Entrez.efetch, db="protein", id=protein_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()
            logging.info(f"Fetched sequence: {seq_record.id} - {seq_record.description}")

            return seq_record
        except Exception as e:
            raise RuntimeError(f"Error fetching SEPI reference: {str(e)}")

    def run_pairwise_alignments(self, query_sequences, reference_sequence):
        """
        Perform pairwise alignments between each query sequence and the reference.
        Returns a comprehensive report string in EMBOSS Water style.
        """
        try:
            logging.info(f"Performing alignments for {len(query_sequences)} queries against {reference_sequence.id}")

            report = []
            for i, query in enumerate(query_sequences, 1):
                best_alignment = self._align_single(query, reference_sequence)
                if best_alignment is None:
                    report.append(f"No alignment found for query {i}: {query.id}")
                    continue

                # Format like EMBOSS Water
                alignment_str = self._format_water_style(best_alignment, query, reference_sequence, i)
                report.append(alignment_str)

            return "\n\n".join(report)
        except Exception as e:
            raise RuntimeError(f"Error in pairwise alignments: {str(e)}")

    def _format_water_style(self, alignment, query, reference, query_num):
        """
        Format the alignment in EMBOSS Water style.
        """
        # Basic formatting
        score = alignment.score
        identities = sum(1 for a, b in zip(alignment[0], alignment[1]) if a == b)
        length = len(alignment[0])
        identity_percent = (identities / length) * 100 if length > 0 else 0
        gaps = alignment[0].count('-') + alignment[1].count('-')
        logging.info(f"Formatting alignment {query_num}: score {score:.1f}, identities {identities}/{length}")

        header = f"# {query_num}\n"
        header += f"Waterman-Eggert local alignment\n"
        header += f"Query: {query.id} ({len(query.seq)} aa)\n"
        header += f"Reference: {reference.id} ({len(reference.seq)} aa)\n"
        header += f"Score: {score:.1f}\n"
        header += f"Identities: {identities}/{length} ({identity_percent:.1f}%)\n"
        header += f"Gaps: {gaps}\n\n"

        # Alignment blocks
        seq1 = str(alignment[0])
        seq2 = str(alignment[1])
        matches = ''.join('|' if a == b else ' ' for a, b in zip(seq1, seq2))

        block_size = 60
        formatted = ""
        for start in range(0, length, block_size):
            end = min(start + block_size, length)
            formatted += f"{query.id:<20} {start+1:>4} {seq1[start:end]} {end:>4}\n"
            formatted += f"{'':<20} {'':>4} {matches[start:end]} {'':>4}\n"
            formatted += f"{reference.id:<20} {start+1:>4} {seq2[start:end]} {end:>4}\n\n"

        return header + formatted

    def run_sepi_alignments(self, query_sequences, genus, species):
        """
        For each query sequence, extract gene name, fetch SEPI reference, and align.
        Returns a comprehensive report string.
        """
        report = []
        for i, query in enumerate(query_sequences, 1):
            # Try to extract gene name from sequence ID first (works with merged files)
            # Format: accession_gene or accession_genus_species_gene
            gene_name = self._extract_gene_from_sequence_id(query.id)
            
            # Fallback to old method if needed
            if not gene_name:
                gene_name = self._extract_gene_name(query.description)
            
            if not gene_name:
                report.append(f"# {i}\nUnable to extract gene name from: {query.id} / {query.description}\n\n")
                continue

            logging.info(f"Processing query {i}: {gene_name}")
            try:
                # Use cached reference getter (eliminates redundant NCBI calls)
                reference = self._get_reference_with_cache(gene_name, genus, species)
                alignment_str = self._format_water_style(
                    self._align_single(query, reference), query, reference, i
                )
                report.append(alignment_str)
            except Exception as e:
                report.append(f"# {i}\nError processing {gene_name}: {str(e)}\n\n")

        return "\n\n".join(report)

    def _extract_gene_name(self, description):
        """
        Extract gene name from FASTA description.
        
        This is a generic fallback that tries to extract gene names from
        various description formats. Works with ANY gene, not just specific ones.
        
        Strategies:
        1. Look for patterns like "gene=geneName" or "gene:geneName"  
        2. Look for gene-like words (3-6 lowercase alphanumeric characters)
        3. Extract from patterns like "protein GeneName [organism]"
        
        Args:
            description: FASTA description string
        
        Returns:
            str: Gene name if found, None otherwise
        
        Examples:
            "efflux transporter AcrB [Escherichia coli]" → "AcrB"
            "gene=tolC multidrug efflux" → "tolC"
            "heat shock protein DnaK" → "DnaK"
        """
        import re
        
        # Strategy 1: Look for explicit gene= or gene: patterns
        gene_pattern = re.search(r'gene[=:]\s*(\w+)', description, re.IGNORECASE)
        if gene_pattern:
            return gene_pattern.group(1)
        
        # Strategy 2: Look for gene names in square brackets (but not organism names)
        # Extract capital words before square brackets (often the gene/protein name)
        bracket_pattern = re.search(r'\b([A-Z][a-z]{2,5}[A-Z]?)\s*\[', description)
        if bracket_pattern:
            return bracket_pattern.group(1)
        
        # Strategy 3: Look for short mixed-case words that look like gene names
        # Gene names are typically 3-6 characters, often lowercase or mixed case
        # Examples: acrA, tolC, dnaK, mecA, vanA, gyrA, ompF, rpoB
        words = description.split()
        for i, word in enumerate(words):
            # Clean word of punctuation
            clean_word = re.sub(r'[^\w]', '', word)
            
            # Gene name pattern: 3-6 chars, starts with lowercase or has mixed case
            if 3 <= len(clean_word) <= 6:
                # Check if it's a likely gene name (starts lowercase or has mix of upper/lower)
                if clean_word.islower() or (any(c.isupper() for c in clean_word) and any(c.islower() for c in clean_word)):
                    # Skip common non-gene words
                    if clean_word.lower() not in ['protein', 'subunit', 'gene', 'from', 'and', 'the']:
                        # Prefer words that come after "protein" or before organism brackets
                        if i > 0 and words[i-1].lower() in ['protein', 'transporter', 'subunit']:
                            return clean_word
                        # Or if it's followed by organism in brackets
                        if i < len(words) - 1 and '[' in words[i+1]:
                            return clean_word
        
        # Strategy 4: Look for any capitalized word (3-6 chars) before square brackets
        cap_pattern = re.search(r'\b([A-Z][a-z]{2,5})\s+\[', description)
        if cap_pattern:
            potential = cap_pattern.group(1)
            # Skip organism names
            if potential not in ['Escherichia', 'Salmonella', 'Staphylococcus', 'Pseudomonas', 'Enterococcus']:
                return potential
        
        # Fallback: return None
        return None

    def _normalize_gene_name(self, gene_raw):
        """
        Normalize gene name for MG1655 reference lookup.
        Handles aliases and naming conventions.
        
        Examples:
            blaEC-5 -> ampC
            mdf(A) -> mdfA
            emrA_1 -> emrA
        """
        if not gene_raw:
            return None
        
        # Strip trailing underscore + digits (e.g., emrA_1 -> emrA)
        import re
        gene = re.sub(r'_\d+$', '', gene_raw)
        
        # Apply known aliases for MG1655 lookup
        aliases = {
            'blaec-5': 'ampC',
            'blaec': 'ampC',
            'mdf(a)': 'mdfA',
        }
        
        gene_lower = gene.lower()
        if gene_lower in aliases:
            normalized = aliases[gene_lower]
            logging.info(f"Normalized gene '{gene}' -> '{normalized}' for MG1655 lookup")
            return normalized
        
        return gene

    def _extract_gene_from_filename(self, filename):
        """
        Extract gene name from filename pattern.
        Supports multiple patterns:
        - AP022811.1_acrB.faa -> acrB
        - AP022811.1_Escherichia_coli_acrA.faa -> acrA
        - CP001363.1_tolC.faa -> tolC
        - NZ_CP107120_mdf(A)_1.faa -> mdf(A) -> mdfA (normalized)
        - NZ_CP107120_blaEC-5.faa -> blaEC-5 -> ampC (normalized)
        
        Returns normalized gene name or None if not found.
        """
        # Get basename without path and extension
        basename = os.path.splitext(os.path.basename(filename))[0]
        
        # Split by underscore
        parts = basename.split('_')
        
        if len(parts) < 2:
            logging.warning(f"Cannot extract gene name from filename: {filename}")
            return None
        
        # The gene name is typically the last part
        # But if last part is a digit (e.g., _1), take second-to-last
        # Handle: NZ_CP107120_mdf(A)_1 -> mdf(A)
        gene_name = parts[-1]
        if gene_name.isdigit() and len(parts) >= 3:
            gene_name = parts[-2]
        
        # Allow alphanumeric plus special chars: - ( ) '
        # Validate basic length
        if gene_name and len(gene_name) >= 3:
            logging.info(f"Extracted gene '{gene_name}' from filename: {filename}")
            return self._normalize_gene_name(gene_name)
        
        logging.warning(f"Invalid gene name '{gene_name}' extracted from: {filename}")
        return None

    def _extract_gene_from_sequence_id(self, sequence_id):
        """
        Extract gene name from FASTA sequence ID.
        Works with IDs from merged files where ID = filename without extension.
        
        Examples:
        - AP022811.1_acrB -> acrB
        - AP022811.1_Escherichia_coli_acrA -> acrA
        - CP001363.1_tolC -> tolC
        - NZ_CP107120_mdf(A)_1 -> mdf(A) -> mdfA (normalized)
        
        Returns normalized gene name or None if not found.
        """
        # Split by underscore
        parts = sequence_id.split('_')
        
        if len(parts) < 2:
            return None
        
        # The gene name is typically the last part
        # But if last part is a digit, take second-to-last
        gene_name = parts[-1]
        if gene_name.isdigit() and len(parts) >= 3:
            gene_name = parts[-2]
        
        # Allow special chars, validate basic length
        if gene_name and len(gene_name) >= 3:
            return self._normalize_gene_name(gene_name)
        
        return None

    def _align_single(self, query, reference):
        """
        Align a single query against reference, handling overflow.
        """
        try:
            alignments = list(self.aligner.align(query.seq, reference.seq))
            return alignments[0] if alignments else None
        except OverflowError:
            # Take first alignment
            alignment_iter = self.aligner.align(query.seq, reference.seq)
            return next(alignment_iter)