"""
FASTA File Analysis - Exercise 3: Basic Sequence Analysis
=========================================================

Learning Objectives:
- Calculate basic sequence statistics (length, composition)
- Work with DNA/RNA nucleotides and amino acids
- Understand GC content and its biological significance
- Practice with dictionaries and mathematical calculations

Background:
Once we can parse FASTA files, we can analyze the sequences. Common analyses include:
- Sequence length
- Nucleotide/amino acid composition
- GC content (for DNA/RNA sequences)
- Finding the most/least common nucleotides or amino acids

GC Content: The percentage of nucleotides that are either Guanine (G) or Cytosine (C).
This is important because GC pairs have stronger hydrogen bonds than AT pairs.

Instructions:
Complete the functions below to analyze DNA sequences from FASTA files.
"""

def calculate_sequence_length(sequence):
    """
    Calculate the length of a sequence.
    
    Args:
        sequence (str): DNA, RNA, or protein sequence
        
    Returns:
        int: Length of the sequence
    """
    # TODO: Return the length of the sequence
    pass


def count_nucleotides(sequence):
    """
    Count the frequency of each nucleotide in a DNA/RNA sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        dict: Dictionary with nucleotide counts
              Keys: 'A', 'T', 'G', 'C' (and 'U' for RNA)
              Values: count of each nucleotide
    """
    # TODO: Count each nucleotide
    # Hint: Use a dictionary to store counts
    # Make sure to handle both uppercase and lowercase
    counts = {}
    # Your code here
    pass


def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA/RNA sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        float: GC content as a percentage (0-100)
    """
    # TODO: Calculate (G + C) / total_nucleotides * 100
    # Handle case where sequence is empty
    pass


def find_most_common_nucleotide(sequence):
    """
    Find the most common nucleotide in a sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        tuple: (nucleotide, count) of the most common nucleotide
    """
    # TODO: Use count_nucleotides and find the maximum
    pass


def analyze_sequence_composition(sequence):
    """
    Perform comprehensive composition analysis of a sequence.
    
    Args:
        sequence (str): DNA, RNA, or protein sequence
        
    Returns:
        dict: Dictionary with analysis results including:
              - length
              - nucleotide_counts (if DNA/RNA)
              - gc_content (if DNA/RNA)
              - most_common
    """
    # TODO: Combine multiple analyses into one comprehensive result
    analysis = {}
    # Your code here
    pass


def compare_gc_content(sequences_dict):
    """
    Compare GC content across multiple sequences.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        
    Returns:
        dict: Dictionary with sequence IDs as keys and GC content as values
    """
    # TODO: Calculate GC content for all sequences
    gc_contents = {}
    # Your code here
    pass


def main():
    """Main function to test your implementations."""
    print("FASTA File Analysis - Exercise 3")
    print("=" * 40)
    
    filename = "data/sample_sequences.fasta"
    
    try:
        # Read and parse the file (using functions from Exercise 2)
        with open(filename, 'r') as file:
            content = file.read()
        
        # Simple parsing function from Exercise 2
        def parse_fasta_simple(fasta_content):
            sequences = []
            lines = fasta_content.strip().split('\n')
            current_header = None
            current_sequence = []
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if current_header is not None:
                        sequences.append((current_header, ''.join(current_sequence)))
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            if current_header is not None:
                sequences.append((current_header, ''.join(current_sequence)))
            
            return sequences
        
        # Parse sequences
        sequences_list = parse_fasta_simple(content)
        
        print(f"Analyzing sequences from: {filename}")
        print()
        
        # Analyze each sequence
        for i, (header, sequence) in enumerate(sequences_list, 1):
            seq_id = header.split()[0][1:]  # Remove '>' and get ID
            print(f"Sequence {i}: {seq_id}")
            print(f"Header: {header}")
            
            # Basic analysis
            length = calculate_sequence_length(sequence)
            print(f"Length: {length}")
            
            # Nucleotide composition
            nucleotide_counts = count_nucleotides(sequence)
            print(f"Nucleotide counts: {nucleotide_counts}")
            
            # GC content
            gc_content = calculate_gc_content(sequence)
            print(f"GC content: {gc_content:.2f}%")
            
            # Most common nucleotide
            most_common = find_most_common_nucleotide(sequence)
            if most_common:
                print(f"Most common nucleotide: {most_common[0]} ({most_common[1]} times)")
            
            # Comprehensive analysis
            analysis = analyze_sequence_composition(sequence)
            print(f"Comprehensive analysis: {analysis}")
            
            print("-" * 50)
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'")
        print("Make sure you're running this script from the correct directory.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()