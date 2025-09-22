import argparse
from Bio import SeqIO
from .engine import AlignmentEngine

def main():
    parser = argparse.ArgumentParser(description="WildTypeAligner: Reference-driven pairwise sequence alignments")
    parser.add_argument('--sequences', required=True, help='Path to multi-FASTA file with query sequences')
    parser.add_argument('--genus', required=True, help='Genus of the organism (e.g., Escherichia)')
    parser.add_argument('--species', required=True, help='Species of the organism (e.g., coli)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--user-reference', help='Path to FASTA file with reference sequence')
    group.add_argument('--sepi', action='store_true', help='Use SEPI to fetch wild-type references for each gene')

    parser.add_argument('--output-file', help='Path to output file (optional, prints to stdout if not provided)')

    args = parser.parse_args()

    # Read query sequences
    try:
        query_sequences = list(SeqIO.parse(args.sequences, 'fasta'))
    except Exception as e:
        print(f"Error reading query sequences: {e}")
        return

    # Run alignments
    engine = AlignmentEngine()
    try:
        if args.user_reference:
            # Single reference for all queries
            reference = SeqIO.read(args.user_reference, 'fasta')
            report = engine.run_pairwise_alignments(query_sequences, reference)
        else:
            # SEPI mode: fetch reference for each query
            report = engine.run_sepi_alignments(query_sequences, args.genus, args.species)
    except Exception as e:
        print(f"Error in alignments: {e}")
        return

    # Output
    if args.output_file:
        try:
            with open(args.output_file, 'w') as f:
                f.write(report)
            print(f"Report saved to {args.output_file}")
        except Exception as e:
            print(f"Error writing to file: {e}")
    else:
        print(report)

if __name__ == '__main__':
    main()