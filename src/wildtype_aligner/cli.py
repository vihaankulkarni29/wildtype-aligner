import argparse
from Bio import SeqIO
from .engine import AlignmentEngine
import os
import glob

def main():
    parser = argparse.ArgumentParser(description="WildTypeAligner: Reference-driven pairwise sequence alignments")
    
    # Input mode: either single file or directory
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--sequences', help='Path to multi-FASTA file with query sequences')
    input_group.add_argument('--input-dir', help='Directory containing .faa files to process')
    
    parser.add_argument('--genus', required=True, help='Genus of the organism (e.g., Escherichia)')
    parser.add_argument('--species', required=True, help='Species of the organism (e.g., coli)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--user-reference', help='Path to FASTA file with reference sequence')
    group.add_argument('--sepi', action='store_true', help='Use SEPI to fetch wild-type references for each gene')

    parser.add_argument('--output-file', help='Path to output file (optional, prints to stdout if not provided)')
    parser.add_argument('--error-log', help='Path to error log file (optional, saves detailed errors if any alignments fail)')
    parser.add_argument('--batch', action='store_true', 
                       help='Enable batch mode: auto-group sequences by gene type and use caching (recommended for large datasets)')
    parser.add_argument('--workers', type=int, default=None,
                       help='Number of parallel worker processes (default: CPU count). Use with --batch for maximum performance')

    # Operational/Entrez configuration
    parser.add_argument('--entrez-email', help='Email for NCBI Entrez (overrides default)')
    parser.add_argument('--ncbi-api-key', help='NCBI API key to increase rate limits')
    parser.add_argument('--entrez-timeout', type=int, default=None, help='Network timeout in seconds for Entrez requests')
    parser.add_argument('--max-retries', type=int, default=3, help='Max retries for Entrez calls on transient failures')
    parser.add_argument('--log-level', choices=['DEBUG','INFO','WARNING','ERROR'], default='INFO', help='Logging verbosity')

    args = parser.parse_args()

    # Read query sequences
    try:
        if args.input_dir:
            # Directory mode: load all .faa files from directory
            if not os.path.isdir(args.input_dir):
                print(f"Error: Directory not found: {args.input_dir}")
                return
            
            # Find all .faa files
            pattern = os.path.join(args.input_dir, "*.faa")
            faa_files = glob.glob(pattern)
            
            if not faa_files:
                print(f"Error: No .faa files found in directory: {args.input_dir}")
                return
            
            print(f"Found {len(faa_files)} .faa files in {args.input_dir}")
            
            # Load all sequences from all files
            query_sequences = []
            for faa_file in sorted(faa_files):
                try:
                    sequences = list(SeqIO.parse(faa_file, 'fasta'))
                    query_sequences.extend(sequences)
                except Exception as e:
                    print(f"Warning: Could not read {faa_file}: {e}")
            
            print(f"Loaded {len(query_sequences)} sequences from {len(faa_files)} files")
        else:
            # Single file mode
            query_sequences = list(SeqIO.parse(args.sequences, 'fasta'))
            print(f"Loaded {len(query_sequences)} sequences from {args.sequences}")
    except Exception as e:
        print(f"Error reading query sequences: {e}")
        return

    # Run alignments
    # Configure logging level early
    import logging
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    # Pass positionally to avoid static analysis issues with keyword params
    engine = AlignmentEngine(
        args.entrez_email,
        args.ncbi_api_key,
        args.entrez_timeout,
        args.max_retries,
    )
    
    # Determine if we should use streaming mode
    # Streaming is beneficial for large datasets to reduce memory usage
    use_streaming = args.output_file is not None and len(query_sequences) >= 50
    
    try:
        if args.user_reference:
            # Single reference for all queries
            reference = SeqIO.read(args.user_reference, 'fasta')
            report = engine.run_pairwise_alignments(query_sequences, reference)
            
            # Write output (no streaming for legacy mode)
            if args.output_file:
                with open(args.output_file, 'w') as f:
                    f.write(report)
                print(f"Report saved to {args.output_file}")
            else:
                print(report)
                
        else:
            # SEPI mode
            if args.batch:
                # Batch mode with optional parallelization and streaming
                if args.workers:
                    # Parallel batch mode (FASTEST!)
                    print(f"Running in PARALLEL BATCH mode with {len(query_sequences)} sequences...")
                    print(f"Using {args.workers} worker processes")
                    
                    if use_streaming:
                        print(f"STREAMING MODE: Writing results directly to {args.output_file}")
                        print("(Memory efficient for large datasets)")
                        with open(args.output_file, 'w') as f:
                            report, summary = engine.run_parallel_batch_alignments(
                                query_sequences, args.genus, args.species, 
                                max_workers=args.workers, output_file=f, error_log=args.error_log
                            )
                        print(f"Report saved to {args.output_file}")
                        print(summary.get_summary())
                    else:
                        report, summary = engine.run_parallel_batch_alignments(
                            query_sequences, args.genus, args.species, 
                            max_workers=args.workers, error_log=args.error_log
                        )
                        if args.output_file:
                            with open(args.output_file, 'w') as f:
                                if report:  # Type guard for Pylance
                                    f.write(report)
                            print(f"Report saved to {args.output_file}")
                        else:
                            if report:
                                print(report)
                        print(summary.get_summary())
                else:
                    # Sequential batch mode (with caching)
                    print(f"Running in BATCH mode with {len(query_sequences)} sequences...")
                    print("Tip: Add --workers flag for parallel processing (faster)")
                    
                    if use_streaming:
                        print(f"STREAMING MODE: Writing results directly to {args.output_file}")
                        with open(args.output_file, 'w') as f:
                            report, summary = engine.run_batch_alignments(
                                query_sequences, args.genus, args.species, 
                                output_file=f, error_log=args.error_log
                            )
                        print(f"Report saved to {args.output_file}")
                        print(summary.get_summary())
                    else:
                        report, summary = engine.run_batch_alignments(
                            query_sequences, args.genus, args.species, error_log=args.error_log
                        )
                        if args.output_file:
                            with open(args.output_file, 'w') as f:
                                if report:
                                    f.write(report)
                            print(f"Report saved to {args.output_file}")
                        else:
                            if report:
                                print(report)
                        print(summary.get_summary())
            else:
                # Legacy mode: fetch reference for each query individually
                print(f"Running in LEGACY mode with {len(query_sequences)} sequences...")
                print("Tip: Use --batch --workers flags for maximum performance")
                report = engine.run_sepi_alignments(query_sequences, args.genus, args.species)
                
                if args.output_file:
                    with open(args.output_file, 'w') as f:
                        f.write(report)
                    print(f"Report saved to {args.output_file}")
                else:
                    print(report)
                    
    except Exception as e:
        print(f"Error in alignments: {e}")
        return

if __name__ == '__main__':
    main()