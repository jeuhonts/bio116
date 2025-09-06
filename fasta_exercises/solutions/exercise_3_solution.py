"""
FASTA File Analysis - Exercise 3 Solution: Basic Sequence Analysis
==================================================================
"""

def calculate_sequence_length(sequence):
    """
    Calculate the length of a sequence.
    
    Args:
        sequence (str): DNA, RNA, or protein sequence
        
    Returns:
        int: Length of the sequence
    """
    return len(sequence)


def count_nucleotides(sequence):
    """
    Count the frequency of each nucleotide in a DNA/RNA sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        dict: Dictionary with nucleotide counts
    """
    # Convert to uppercase for consistent counting
    sequence = sequence.upper()
    
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'U': 0}
    
    for nucleotide in sequence:
        if nucleotide in counts:
            counts[nucleotide] += 1
    
    return counts


def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA/RNA sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        float: GC content as a percentage (0-100)
    """
    if len(sequence) == 0:
        return 0.0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_nucleotides = len(sequence)
    
    return (gc_count / total_nucleotides) * 100


def find_most_common_nucleotide(sequence):
    """
    Find the most common nucleotide in a sequence.
    
    Args:
        sequence (str): DNA or RNA sequence
        
    Returns:
        tuple: (nucleotide, count) of the most common nucleotide
    """
    counts = count_nucleotides(sequence)
    
    # Find the nucleotide with the maximum count
    max_nucleotide = max(counts, key=counts.get)
    max_count = counts[max_nucleotide]
    
    return (max_nucleotide, max_count)


def analyze_sequence_composition(sequence):
    """
    Perform comprehensive composition analysis of a sequence.
    
    Args:
        sequence (str): DNA, RNA, or protein sequence
        
    Returns:
        dict: Dictionary with analysis results
    """
    analysis = {
        'length': calculate_sequence_length(sequence),
        'nucleotide_counts': count_nucleotides(sequence),
        'gc_content': calculate_gc_content(sequence),
        'most_common': find_most_common_nucleotide(sequence)
    }
    
    return analysis


def compare_gc_content(sequences_dict):
    """
    Compare GC content across multiple sequences.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        
    Returns:
        dict: Dictionary with sequence IDs as keys and GC content as values
    """
    gc_contents = {}
    
    for seq_id, seq_info in sequences_dict.items():
        sequence = seq_info['sequence']
        gc_contents[seq_id] = calculate_gc_content(sequence)
    
    return gc_contents


def main():
    """Main function to test the implementations."""
    print("FASTA File Analysis - Exercise 3 Solution")
    print("=" * 45)
    
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
            print(f"Comprehensive analysis keys: {list(analysis.keys())}")
            
            print("-" * 50)
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'")
        print("Make sure you're running this script from the correct directory.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()