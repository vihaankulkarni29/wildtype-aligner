# WildTypeAligner

## What is WildTypeAligner?

WildTypeAligner is a specialized bioinformatics tool designed for researchers studying antimicrobial resistance. It helps you compare protein sequences from resistant bacteria against their "normal" (wild-type) counterparts to identify mutations that cause drug resistance.

**Imagine you have bacterial proteins that have developed resistance to antibiotics. This tool automatically finds the original, non-resistant versions of those same proteins from the same bacterial species and shows you exactly where the differences are.**

## Key Features

- 🔬 **Automatic Reference Finding**: No need to manually search for reference sequences - the tool fetches them directly from NCBI databases
- 🎯 **Species-Specific Matching**: Ensures comparisons are made against proteins from the exact same bacterial species
- 📊 **Clear Mutation Reports**: Shows alignments in an easy-to-read format with scores, identities, and gap information
- ⚡ **Batch Processing**: Analyze multiple proteins at once
- 🧬 **Works with Any Bacteria**: From E. coli to Mycobacterium tuberculosis - any bacterial species

## Who Should Use This Tool?

- **Microbiologists** studying antibiotic resistance
- **Infectious Disease Researchers** tracking bacterial evolution
- **Clinical Laboratories** identifying resistance mechanisms
- **Bioinformatics Students** learning sequence analysis

## Prerequisites

Before using WildTypeAligner, you need:

1. **Python 3.8 or higher** installed on your computer
2. **Internet connection** (to fetch reference sequences from NCBI)
3. **Basic knowledge of FASTA file format** (standard for biological sequences)

## Installation

### Step 1: Download the Tool

```bash
# Clone the repository from GitHub
git clone https://github.com/vihaankulkarni29/wildtype-aligner.git
cd wildtype-aligner
```

### Step 2: Install the Software

```bash
# Install the tool and its dependencies
pip install .
```

**That's it!** The tool is now ready to use.

## How to Use WildTypeAligner

### What You Need to Prepare

1. **A FASTA file** containing your protein sequences (the resistant ones you want to analyze)
2. **Know your bacterial species** (genus and species name, like "Escherichia coli")

### Two Ways to Use the Tool

#### Method 1: Automatic Mode (Recommended)

Let the tool automatically find reference sequences:

```bash
wildtype-aligner --sequences your_proteins.fasta --genus Escherichia --species coli --sepi --output-file results.txt
```

#### Method 2: Manual Reference Mode

If you have your own reference sequences:

```bash
wildtype-aligner --sequences your_proteins.fasta --genus Escherichia --species coli --user-reference reference_protein.fasta --output-file results.txt
```

### Understanding the Results

The output file contains alignment reports like this:

```
# 1
Waterman-Eggert local alignment
Query: resistant_protein_1 (400 aa)
Reference: wild_type_protein (395 aa)
Score: 350.0
Identities: 350/400 (87.5%)
Gaps: 5

Query_Sequence    1 MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKT   60
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||
Reference_Seq     1 MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKT   60

... (continues showing differences)
```

**What the numbers mean:**
- **Score**: How well the sequences match (higher = better match)
- **Identities**: Percentage of amino acids that are identical
- **Gaps**: Insertions/deletions in the sequence
- **The alignment**: Shows exactly where mutations occur (differences marked)

### Example Workflow

Let's say you have E. coli proteins that are resistant to tetracycline:

1. **Prepare your data**: Create `resistant_proteins.fasta` with your resistant protein sequences
2. **Run the analysis**:
   ```bash
   wildtype-aligner --sequences resistant_proteins.fasta --genus Escherichia --species coli --sepi --output-file tetracycline_resistance.txt
   ```
3. **Interpret results**: Look for mutations in genes like tetA, tetR, etc.

## Troubleshooting

### "No protein sequence found"
- Check that genus and species names are spelled correctly
- Ensure gene names in your FASTA headers are standard (e.g., "AcrB" not "acrB_v2")

### "Internet connection error"
- Make sure you have internet access for NCBI queries
- NCBI has usage limits; wait a few minutes if you get rate-limited

### "FASTA file format error"
- Ensure your input file follows FASTA format:
  ```
  >protein_name description
  MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKT...
  ```

## Scientific Background

Antimicrobial resistance occurs when bacteria evolve to survive antibiotic treatment. This often happens through mutations in proteins that:
- Pump out antibiotics (efflux pumps like AcrAB-TolC)
- Modify antibiotics (enzymes)
- Change antibiotic targets (ribosomes, cell walls)

WildTypeAligner helps identify these mutations by comparing resistant proteins against their wild-type ancestors, revealing the genetic changes responsible for resistance.

## Citation

If you use WildTypeAligner in your research, please cite:

```
Kulkarni V. WildTypeAligner: Reference-driven pairwise sequence alignment tool for antimicrobial resistance research. 2025.
```

## Support

For questions or issues:
- Check the [Issues](https://github.com/vihaankulkarni29/wildtype-aligner/issues) page on GitHub
- Email: vihaan.kulkarni@wildtypealigner.com

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.