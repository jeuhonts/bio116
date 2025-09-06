# BIO116 - Bioinformatics Course

This repository contains Python coding exercises for learning bioinformatics, specifically focused on FASTA sequence file analysis.

## Course Content

### FASTA File Analysis Exercises

A comprehensive series of Python exercises designed to teach students how to work with biological sequence data in FASTA format. These exercises progress from basic file operations to advanced sequence analysis techniques.

#### What You'll Learn

- Understanding FASTA file format structure
- Reading and parsing biological sequence files
- Basic sequence analysis (length, composition, GC content)
- Advanced sequence analysis (translation, ORF finding, reverse complement)
- Multi-sequence comparison and analysis
- Bioinformatics programming best practices

#### Exercise Structure

Each exercise includes:
- **Learning objectives** - Clear goals for what you'll accomplish
- **Background information** - Biological and computational context
- **Starter code** - Functions to complete with TODO comments
- **Complete solutions** - Reference implementations for learning
- **Test data** - Sample FASTA files for practice

## Getting Started

### Prerequisites

- Python 3.6 or higher
- Basic understanding of Python programming (variables, functions, loops)
- Text editor or IDE for Python development

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/jeuhonts/bio116.git
   cd bio116
   ```

2. Navigate to the exercises directory:
   ```bash
   cd fasta_exercises
   ```

### Running the Exercises

Each exercise can be run independently:

```bash
# Work on Exercise 1 - Basic file reading
python3 exercise_1_basic_reading.py

# Check your work against the solution
python3 solutions/exercise_1_solution.py
```

## Exercise Overview

### Exercise 1: Basic File Reading and FASTA Format
**File:** `exercise_1_basic_reading.py`

Learn the fundamentals of:
- Opening and reading files in Python
- Understanding FASTA format structure
- Counting sequences and extracting headers
- Basic string manipulation

**Key Functions:**
- `read_fasta_file()` - Read file contents
- `count_sequences()` - Count FASTA entries
- `get_sequence_headers()` - Extract sequence headers

### Exercise 2: Parsing Headers and Sequences
**File:** `exercise_2_parsing.py`

Build on basic file reading to:
- Parse FASTA into structured data
- Separate headers from sequences
- Handle multi-line sequences
- Create dictionaries for data storage

**Key Functions:**
- `parse_fasta_simple()` - Basic FASTA parsing
- `parse_fasta_dict()` - Structured data parsing
- `get_sequence_by_id()` - Retrieve specific sequences

### Exercise 3: Basic Sequence Analysis
**File:** `exercise_3_analysis.py`

Perform fundamental sequence analysis:
- Calculate sequence statistics
- Analyze nucleotide composition
- Calculate GC content
- Find most common nucleotides

**Key Functions:**
- `calculate_sequence_length()` - Get sequence length
- `count_nucleotides()` - Count A, T, G, C
- `calculate_gc_content()` - Calculate GC percentage
- `analyze_sequence_composition()` - Comprehensive analysis

### Exercise 4: Advanced Sequence Analysis
**File:** `exercise_4_advanced.py`

Advanced bioinformatics techniques:
- Pattern finding in sequences
- DNA to protein translation
- Open Reading Frame (ORF) identification
- Reverse complement calculation
- Molecular weight estimation

**Key Functions:**
- `find_pattern()` - Find sequence motifs
- `translate_dna()` - Genetic code translation
- `find_orfs()` - Identify protein-coding regions
- `reverse_complement()` - Calculate reverse complement

### Exercise 5: Working with Multiple Sequences
**File:** `exercise_5_multiple.py`

Multi-sequence analysis:
- Compare multiple sequences
- Calculate dataset statistics
- Find common motifs across sequences
- Create similarity matrices
- Generate comprehensive reports

**Key Functions:**
- `compare_sequence_lengths()` - Length statistics
- `find_common_motifs()` - Cross-sequence motif analysis
- `create_similarity_matrix()` - Sequence comparison
- `generate_sequence_report()` - Comprehensive reporting

## Sample Data

The exercises include sample FASTA files:

- **`data/sample_sequences.fasta`** - DNA sequences from various organisms
  - Human insulin gene
  - E. coli beta-galactosidase fragment
  - Yeast alcohol dehydrogenase
  - Arabidopsis chlorophyll binding protein

- **`data/protein_sequences.fasta`** - Protein sequences
  - Human hemoglobin alpha chain
  - E. coli Trp operon peptide
  - Yeast cytochrome c

## Tips for Success

1. **Start Simple**: Begin with Exercise 1 and work through them sequentially
2. **Read the Background**: Each exercise includes important biological context
3. **Test Frequently**: Run your code often to catch errors early
4. **Use the Solutions**: Compare your work with provided solutions
5. **Experiment**: Try the functions with different data and parameters

## Common Pitfalls

- **Case Sensitivity**: DNA sequences can be uppercase or lowercase
- **File Paths**: Make sure you're in the correct directory when running scripts
- **Empty Sequences**: Handle edge cases like empty files or sequences
- **Multi-line Sequences**: FASTA sequences often span multiple lines

## Additional Resources

- [FASTA Format Specification](https://en.wikipedia.org/wiki/FASTA_format)
- [Genetic Code](https://en.wikipedia.org/wiki/Genetic_code)
- [BioPython Documentation](https://biopython.org/) - For more advanced work
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) - For sequence searching

## Contributing

This is an educational repository. If you find bugs or have suggestions for improvements:

1. Check existing issues first
2. Create detailed bug reports
3. Suggest improvements to exercises or documentation

## License

This educational content is provided for learning purposes. Please respect any licensing terms for sequence data from external sources.
