"""
FASTA File Analysis - Exercise 5: Working with Multiple Sequences
=================================================================

Learning Objectives:
- Compare multiple sequences
- Calculate statistics across datasets
- Find consensus sequences
- Perform sequence alignments (basic)
- Generate reports and summaries

Background:
When working with multiple sequences, we often want to compare them, find similarities
and differences, calculate statistics across the entire dataset, and generate
comprehensive reports. This exercise focuses on multi-sequence analysis.

Instructions:
Complete the functions below to analyze multiple sequences from FASTA files.
"""

def load_sequences_from_file(filename):
    """
    Load all sequences from a FASTA file into a dictionary.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        dict: Dictionary with sequence IDs as keys and sequence info as values
    """
    # TODO: Load and parse the FASTA file
    # Use parsing functions from previous exercises
    sequences = {}
    # Your code here
    pass


def compare_sequence_lengths(sequences_dict):
    """
    Compare lengths of all sequences in the dataset.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: Statistics including min, max, average length
    """
    # TODO: Calculate length statistics
    stats = {
        'min_length': 0,
        'max_length': 0,
        'average_length': 0,
        'total_sequences': 0
    }
    # Your code here
    pass


def calculate_gc_content_distribution(sequences_dict):
    """
    Calculate GC content for all sequences and return distribution statistics.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: GC content statistics and individual values
    """
    # TODO: Calculate GC content for each sequence and overall statistics
    result = {
        'individual_gc': {},  # seq_id: gc_content
        'min_gc': 0,
        'max_gc': 0,
        'average_gc': 0
    }
    # Your code here
    pass


def find_common_motifs(sequences_dict, motif_length=3):
    """
    Find motifs that appear in multiple sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        motif_length (int): Length of motifs to search for
        
    Returns:
        dict: Dictionary with motifs as keys and lists of sequence IDs as values
    """
    # TODO: Find motifs that appear in multiple sequences
    # Look for all possible motifs of given length
    # Count how many sequences contain each motif
    common_motifs = {}
    # Your code here
    pass


def calculate_sequence_similarity(seq1, seq2):
    """
    Calculate similarity between two sequences as percentage of identical positions.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        float: Similarity percentage (0-100)
    """
    # TODO: Calculate similarity
    # Compare sequences position by position
    # Handle sequences of different lengths
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0
    
    # Your code here
    pass


def create_similarity_matrix(sequences_dict):
    """
    Create a similarity matrix comparing all sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: Nested dictionary with similarity values
    """
    # TODO: Calculate similarity between all pairs of sequences
    similarity_matrix = {}
    # Your code here
    pass


def find_consensus_nucleotides(sequences_dict, position):
    """
    Find the most common nucleotide at a specific position across all sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        position (int): Position to analyze (0-based)
        
    Returns:
        tuple: (most_common_nucleotide, frequency)
    """
    # TODO: Find consensus at specific position
    # Count nucleotides at the given position across all sequences
    nucleotide_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    # Your code here
    pass


def generate_sequence_report(sequences_dict):
    """
    Generate a comprehensive report about the sequence dataset.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        str: Formatted report string
    """
    # TODO: Create a comprehensive report
    # Include length statistics, GC content, composition, etc.
    report = "FASTA Dataset Analysis Report\n"
    report += "=" * 40 + "\n\n"
    
    # Your code here - add various statistics and analyses
    
    return report


def main():
    """Main function to test your implementations."""
    print("FASTA File Analysis - Exercise 5")
    print("=" * 40)
    
    # Test with both DNA and protein sequences
    filenames = ["data/sample_sequences.fasta", "data/protein_sequences.fasta"]
    
    for filename in filenames:
        try:
            print(f"\nAnalyzing file: {filename}")
            print("-" * 30)
            
            # Load sequences
            sequences = load_sequences_from_file(filename)
            print(f"Loaded {len(sequences)} sequences")
            
            # Compare lengths
            length_stats = compare_sequence_lengths(sequences)
            print(f"Length statistics: {length_stats}")
            
            # GC content analysis (for DNA sequences)
            if "sample_sequences" in filename:  # DNA sequences
                gc_stats = calculate_gc_content_distribution(sequences)
                print(f"GC content range: {gc_stats['min_gc']:.2f}% - {gc_stats['max_gc']:.2f}%")
                print(f"Average GC content: {gc_stats['average_gc']:.2f}%")
            
            # Find common motifs
            common_motifs = find_common_motifs(sequences, motif_length=3)
            print(f"Found {len(common_motifs)} common motifs (length 3)")
            
            # Show some common motifs
            motif_items = list(common_motifs.items())
            for motif, seq_ids in motif_items[:5]:  # Show first 5
                print(f"  {motif}: appears in {len(seq_ids)} sequences")
            
            # Similarity matrix
            similarity_matrix = create_similarity_matrix(sequences)
            print("Similarity matrix (first few comparisons):")
            seq_ids = list(sequences.keys())
            for i, seq1 in enumerate(seq_ids[:2]):
                for j, seq2 in enumerate(seq_ids[:2]):
                    if seq1 in similarity_matrix and seq2 in similarity_matrix[seq1]:
                        sim = similarity_matrix[seq1][seq2]
                        print(f"  {seq1} vs {seq2}: {sim:.2f}%")
            
            # Generate comprehensive report
            report = generate_sequence_report(sequences)
            print("\nComprehensive Report:")
            print(report[:500] + "..." if len(report) > 500 else report)
            
        except FileNotFoundError:
            print(f"Warning: Could not find file '{filename}' - skipping")
        except Exception as e:
            print(f"Error analyzing {filename}: {e}")


if __name__ == "__main__":
    main()