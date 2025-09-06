# FASTA Exercises - Quick Start Guide

This guide helps you get started with the FASTA analysis exercises quickly.

## Quick Setup

1. Navigate to the exercises directory:
   ```bash
   cd fasta_exercises
   ```

2. List available exercises:
   ```bash
   ls exercise_*.py
   ```

3. Start with Exercise 1:
   ```bash
   python3 exercise_1_basic_reading.py
   ```

## Exercise Progression

Work through the exercises in order for the best learning experience:

1. **Exercise 1** → Basic file operations and FASTA format
2. **Exercise 2** → Parsing and data structures  
3. **Exercise 3** → Sequence analysis and statistics
4. **Exercise 4** → Advanced bioinformatics techniques
5. **Exercise 5** → Multi-sequence analysis

## File Structure

```
fasta_exercises/
├── data/                          # Sample FASTA files
│   ├── sample_sequences.fasta     # DNA sequences for practice
│   └── protein_sequences.fasta    # Protein sequences for practice
├── solutions/                     # Complete solutions
│   ├── exercise_1_solution.py
│   ├── exercise_2_solution.py
│   ├── exercise_3_solution.py
│   ├── exercise_4_solution.py
│   └── exercise_5_solution.py
├── exercise_1_basic_reading.py    # Basic file reading
├── exercise_2_parsing.py          # FASTA parsing
├── exercise_3_analysis.py         # Sequence analysis
├── exercise_4_advanced.py         # Advanced analysis
└── exercise_5_multiple.py         # Multi-sequence work
```

## How to Work on Exercises

### Step 1: Read the Exercise
Each exercise starts with:
- Learning objectives
- Background information
- Instructions

### Step 2: Complete the Functions
Look for `# TODO:` comments and implement the missing code:

```python
def count_nucleotides(sequence):
    """Count frequency of each nucleotide."""
    # TODO: Count each nucleotide
    counts = {}
    # Your code here
    pass
```

### Step 3: Test Your Code
Run the exercise to see if your implementation works:

```bash
python3 exercise_1_basic_reading.py
```

### Step 4: Compare with Solutions
Check your work against the provided solutions:

```bash
python3 solutions/exercise_1_solution.py
```

## Key Concepts

### FASTA Format
```
>sequence_1 Human insulin gene
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTG
GGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCAC
```

- **Header line**: Starts with `>`, contains ID and description
- **Sequence lines**: Actual sequence data (can span multiple lines)

### Important Python Concepts

**File Reading:**
```python
with open(filename, 'r') as file:
    content = file.read()
```

**String Methods:**
```python
# Split into lines
lines = content.split('\n')

# Check if line starts with '>'
if line.startswith('>'):
    # This is a header line

# Convert to uppercase
sequence = sequence.upper()

# Count specific characters
gc_count = sequence.count('G') + sequence.count('C')
```

**Dictionaries:**
```python
# Store sequence information
sequences = {
    'seq_1': {
        'description': 'Human gene',
        'sequence': 'ATGGCC...'
    }
}
```

## Common Functions You'll Implement

### Exercise 1 - Basic Operations
- Read files
- Count sequences
- Extract headers

### Exercise 2 - Data Parsing
- Parse FASTA format
- Create data structures
- Retrieve sequences by ID

### Exercise 3 - Sequence Analysis
- Calculate length and composition
- Compute GC content
- Find most common nucleotides

### Exercise 4 - Advanced Analysis
- Find patterns and motifs
- Translate DNA to protein
- Find Open Reading Frames
- Calculate reverse complement

### Exercise 5 - Multi-Sequence Analysis
- Compare multiple sequences
- Find common motifs
- Calculate similarity
- Generate reports

## Tips for Success

1. **Read Error Messages**: Python error messages are helpful
2. **Print Variables**: Use `print()` to debug your code
3. **Start Small**: Test with simple examples first
4. **Use the Python REPL**: Try small code snippets interactively

## Getting Help

If you get stuck:

1. **Re-read the instructions** and background information
2. **Check the sample data** in the `data/` directory
3. **Look at solutions** for similar functions in earlier exercises
4. **Test with simple inputs** to understand the problem
5. **Use Python's help system**: `help(function_name)`

## Sample Session

Here's what a typical working session might look like:

```bash
# Navigate to exercises
cd fasta_exercises

# Work on Exercise 1
python3 exercise_1_basic_reading.py
# Error: FileNotFoundError

# Check if data file exists
ls data/
# sample_sequences.fasta  protein_sequences.fasta

# Run again - now it works!
python3 exercise_1_basic_reading.py
# FASTA File Analysis - Exercise 1
# ================================
# Successfully read file: data/sample_sequences.fasta
# ...

# Check solution
python3 solutions/exercise_1_solution.py
# Compare output with your implementation

# Move to next exercise
python3 exercise_2_parsing.py
```

## Next Steps

After completing all exercises, you might want to:

1. **Modify the sample data** - Try your own FASTA files
2. **Extend the functions** - Add new features or analyses
3. **Combine exercises** - Create a comprehensive analysis script
4. **Explore BioPython** - A powerful library for bioinformatics
5. **Learn about sequence databases** - NCBI, UniProt, etc.

Happy coding!