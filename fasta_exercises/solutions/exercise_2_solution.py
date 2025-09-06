"""
FASTA File Analysis - Exercise 2 Solution: Parsing Headers and Sequences
=========================================================================
"""

def parse_fasta_simple(fasta_content):
    """
    Parse FASTA content into a list of (header, sequence) tuples.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        list: List of tuples where each tuple is (header, sequence)
              header includes the '>' character
              sequence is concatenated from all sequence lines
    """
    sequences = []
    lines = fasta_content.strip().split('\n')
    
    current_header = None
    current_sequence = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence if it exists
            if current_header is not None:
                sequences.append((current_header, ''.join(current_sequence)))
            
            # Start new sequence
            current_header = line
            current_sequence = []
        else:
            # Add to current sequence
            current_sequence.append(line)
    
    # Don't forget the last sequence
    if current_header is not None:
        sequences.append((current_header, ''.join(current_sequence)))
    
    return sequences


def parse_fasta_dict(fasta_content):
    """
    Parse FASTA content into a dictionary with sequence IDs as keys.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        dict: Dictionary where keys are sequence IDs (without '>') 
              and values are dictionaries with 'description' and 'sequence' keys
    """
    result = {}
    sequences_list = parse_fasta_simple(fasta_content)
    
    for header, sequence in sequences_list:
        # Remove '>' and split header into ID and description
        header_without_gt = header[1:]  # Remove '>'
        parts = header_without_gt.split(' ', 1)  # Split on first space
        
        seq_id = parts[0]
        description = parts[1] if len(parts) > 1 else ""
        
        result[seq_id] = {
            'description': description,
            'sequence': sequence
        }
    
    return result


def get_sequence_by_id(sequences_dict, sequence_id):
    """
    Retrieve a specific sequence by its ID.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        sequence_id (str): ID of the sequence to retrieve
        
    Returns:
        dict: Dictionary with sequence information, or None if not found
    """
    return sequences_dict.get(sequence_id)


def get_all_sequence_ids(sequences_dict):
    """
    Get a list of all sequence IDs in the dataset.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        
    Returns:
        list: List of all sequence IDs
    """
    return list(sequences_dict.keys())


def main():
    """Main function to test the implementations."""
    print("FASTA File Analysis - Exercise 2 Solution")
    print("=" * 45)
    
    filename = "data/sample_sequences.fasta"
    
    try:
        # Read the file
        with open(filename, 'r') as file:
            content = file.read()
        
        print(f"Parsing file: {filename}")
        print()
        
        # Test simple parsing
        sequences_list = parse_fasta_simple(content)
        print("Simple parsing results:")
        print(f"Found {len(sequences_list)} sequences")
        for i, (header, sequence) in enumerate(sequences_list, 1):
            print(f"  {i}. {header}")
            print(f"     Sequence length: {len(sequence)}")
            print(f"     First 50 chars: {sequence[:50]}...")
        print()
        
        # Test dictionary parsing
        sequences_dict = parse_fasta_dict(content)
        print("Dictionary parsing results:")
        print(f"Found {len(sequences_dict)} sequences")
        
        # Get all IDs
        all_ids = get_all_sequence_ids(sequences_dict)
        print(f"Sequence IDs: {all_ids}")
        print()
        
        # Test retrieving specific sequences
        if all_ids:
            first_id = all_ids[0]
            seq_info = get_sequence_by_id(sequences_dict, first_id)
            if seq_info:
                print(f"Details for '{first_id}':")
                print(f"  Description: {seq_info['description']}")
                print(f"  Sequence length: {len(seq_info['sequence'])}")
                print(f"  First 50 chars: {seq_info['sequence'][:50]}...")
            
        print()
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'")
        print("Make sure you're running this script from the correct directory.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()