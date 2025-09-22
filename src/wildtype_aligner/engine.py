from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
import textwrap
import logging

logging.basicConfig(level=logging.INFO)

class AlignmentEngine:
    def __init__(self):
        # Set email for NCBI Entrez
        Entrez.email = "vihaan.kulkarni@wildtypealigner.com"  # Replace with your actual email
        # Initialize aligner
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # Use global alignment
        logging.info("AlignmentEngine initialized")

    def _fetch_sepi_reference(self, gene_name, genus, species):
        """
        Fetch the canonical protein sequence from NCBI for the given gene, genus, and species.
        Returns a SeqRecord object.
        """
        try:
            # Construct search term with genus and species for precision
            search_term = f"{gene_name}[Gene] AND {genus} {species}[Organism]"
            logging.info(f"Searching NCBI for {gene_name} in {genus} {species}")
            # Search for protein records
            handle = Entrez.esearch(db="protein", term=search_term, retmax=5)  # Get a few to choose the best
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                raise ValueError(f"No protein sequence found for gene {gene_name} in {genus} {species}")

            # Take the first result (assuming it's the canonical one)
            protein_id = record["IdList"][0]
            # Fetch the sequence
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
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
            # Extract gene name from description
            gene_name = self._extract_gene_name(query.description)
            if not gene_name:
                report.append(f"# {i}\nUnable to extract gene name from: {query.description}\n\n")
                continue

            logging.info(f"Processing query {i}: {gene_name}")
            try:
                reference = self._fetch_sepi_reference(gene_name, genus, species)
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
        Assumes format like "... subunit AcrA ..."
        """
        # Simple extraction: look for common gene patterns
        words = description.split()
        for word in words:
            if word.startswith('Acr') and len(word) > 3:  # AcrA, AcrB, etc.
                return word
        # Fallback: return None
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