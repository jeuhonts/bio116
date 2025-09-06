"""
FASTA File Analysis - Exercise 5 Solution: Working with Multiple Sequences
==========================================================================
"""

def load_sequences_from_file(filename):
    """
    Load all sequences from a FASTA file into a dictionary.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        dict: Dictionary with sequence IDs as keys and sequence info as values
    """
    with open(filename, 'r') as file:
        content = file.read()
    
    # Parse FASTA content
    sequences = {}
    lines = content.strip().split('\n')
    current_header = None
    current_sequence = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence
            if current_header is not None:
                header_parts = current_header[1:].split(' ', 1)
                seq_id = header_parts[0]
                description = header_parts[1] if len(header_parts) > 1 else ""
                sequences[seq_id] = {
                    'description': description,
                    'sequence': ''.join(current_sequence)
                }
            
            # Start new sequence
            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)
    
    # Don't forget the last sequence
    if current_header is not None:
        header_parts = current_header[1:].split(' ', 1)
        seq_id = header_parts[0]
        description = header_parts[1] if len(header_parts) > 1 else ""
        sequences[seq_id] = {
            'description': description,
            'sequence': ''.join(current_sequence)
        }
    
    return sequences


def compare_sequence_lengths(sequences_dict):
    """
    Compare lengths of all sequences in the dataset.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: Statistics including min, max, average length
    """
    if not sequences_dict:
        return {'min_length': 0, 'max_length': 0, 'average_length': 0, 'total_sequences': 0}
    
    lengths = [len(seq_info['sequence']) for seq_info in sequences_dict.values()]
    
    stats = {
        'min_length': min(lengths),
        'max_length': max(lengths),
        'average_length': sum(lengths) / len(lengths),
        'total_sequences': len(sequences_dict)
    }
    
    return stats


def calculate_gc_content_distribution(sequences_dict):
    """
    Calculate GC content for all sequences and return distribution statistics.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: GC content statistics and individual values
    """
    def calculate_gc_content(sequence):
        if len(sequence) == 0:
            return 0.0
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    individual_gc = {}
    gc_values = []
    
    for seq_id, seq_info in sequences_dict.items():
        gc_content = calculate_gc_content(seq_info['sequence'])
        individual_gc[seq_id] = gc_content
        gc_values.append(gc_content)
    
    result = {
        'individual_gc': individual_gc,
        'min_gc': min(gc_values) if gc_values else 0,
        'max_gc': max(gc_values) if gc_values else 0,
        'average_gc': sum(gc_values) / len(gc_values) if gc_values else 0
    }
    
    return result


def find_common_motifs(sequences_dict, motif_length=3):
    """
    Find motifs that appear in multiple sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        motif_length (int): Length of motifs to search for
        
    Returns:
        dict: Dictionary with motifs as keys and lists of sequence IDs as values
    """
    motif_occurrences = {}
    
    for seq_id, seq_info in sequences_dict.items():
        sequence = seq_info['sequence'].upper()
        
        # Find all motifs of given length in this sequence
        sequence_motifs = set()
        for i in range(len(sequence) - motif_length + 1):
            motif = sequence[i:i+motif_length]
            sequence_motifs.add(motif)
        
        # Add this sequence to the motif occurrence lists
        for motif in sequence_motifs:
            if motif not in motif_occurrences:
                motif_occurrences[motif] = []
            motif_occurrences[motif].append(seq_id)
    
    # Filter to only motifs that appear in multiple sequences
    common_motifs = {motif: seq_ids for motif, seq_ids in motif_occurrences.items() 
                    if len(seq_ids) > 1}
    
    return common_motifs


def calculate_sequence_similarity(seq1, seq2):
    """
    Calculate similarity between two sequences as percentage of identical positions.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        float: Similarity percentage (0-100)
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0
    
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Compare up to the length of the shorter sequence
    min_length = min(len(seq1), len(seq2))
    matches = 0
    
    for i in range(min_length):
        if seq1[i] == seq2[i]:
            matches += 1
    
    return (matches / min_length) * 100


def create_similarity_matrix(sequences_dict):
    """
    Create a similarity matrix comparing all sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        dict: Nested dictionary with similarity values
    """
    similarity_matrix = {}
    seq_ids = list(sequences_dict.keys())
    
    for seq1_id in seq_ids:
        similarity_matrix[seq1_id] = {}
        for seq2_id in seq_ids:
            seq1 = sequences_dict[seq1_id]['sequence']
            seq2 = sequences_dict[seq2_id]['sequence']
            similarity = calculate_sequence_similarity(seq1, seq2)
            similarity_matrix[seq1_id][seq2_id] = similarity
    
    return similarity_matrix


def find_consensus_nucleotides(sequences_dict, position):
    """
    Find the most common nucleotide at a specific position across all sequences.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        position (int): Position to analyze (0-based)
        
    Returns:
        tuple: (most_common_nucleotide, frequency)
    """
    nucleotide_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    total_sequences = 0
    
    for seq_info in sequences_dict.values():
        sequence = seq_info['sequence'].upper()
        if position < len(sequence):
            nucleotide = sequence[position]
            if nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] += 1
                total_sequences += 1
    
    if total_sequences == 0:
        return ('N', 0)
    
    # Find most common nucleotide
    most_common = max(nucleotide_counts, key=nucleotide_counts.get)
    frequency = nucleotide_counts[most_common] / total_sequences
    
    return (most_common, frequency)


def generate_sequence_report(sequences_dict):
    """
    Generate a comprehensive report about the sequence dataset.
    
    Args:
        sequences_dict (dict): Dictionary of sequences
        
    Returns:
        str: Formatted report string
    """
    report = "FASTA Dataset Analysis Report\n"
    report += "=" * 40 + "\n\n"
    
    if not sequences_dict:
        report += "No sequences found in dataset.\n"
        return report
    
    # Basic statistics
    length_stats = compare_sequence_lengths(sequences_dict)
    report += f"Dataset Overview:\n"
    report += f"  Total sequences: {length_stats['total_sequences']}\n"
    report += f"  Length range: {length_stats['min_length']} - {length_stats['max_length']} nt\n"
    report += f"  Average length: {length_stats['average_length']:.2f} nt\n\n"
    
    # Individual sequence details
    report += "Individual Sequences:\n"
    for seq_id, seq_info in sequences_dict.items():
        sequence = seq_info['sequence']
        report += f"  {seq_id}: {len(sequence)} nt - {seq_info['description']}\n"
    
    # GC content analysis (assuming DNA sequences)
    try:
        gc_stats = calculate_gc_content_distribution(sequences_dict)
        report += f"\nGC Content Analysis:\n"
        report += f"  Range: {gc_stats['min_gc']:.2f}% - {gc_stats['max_gc']:.2f}%\n"
        report += f"  Average: {gc_stats['average_gc']:.2f}%\n"
        
        for seq_id, gc_content in gc_stats['individual_gc'].items():
            report += f"    {seq_id}: {gc_content:.2f}%\n"
    except:
        pass  # Skip if not DNA sequences
    
    # Common motifs
    common_motifs = find_common_motifs(sequences_dict, motif_length=3)
    report += f"\nCommon Motifs (length 3): {len(common_motifs)} found\n"
    
    # Show top 5 most common motifs
    motif_items = sorted(common_motifs.items(), key=lambda x: len(x[1]), reverse=True)
    for motif, seq_ids in motif_items[:5]:
        report += f"  {motif}: appears in {len(seq_ids)} sequences\n"
    
    return report


def main():
    """Main function to test the implementations."""
    print("FASTA File Analysis - Exercise 5 Solution")
    print("=" * 45)
    
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