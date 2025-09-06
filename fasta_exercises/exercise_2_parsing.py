"""
FASTA File Analysis - Exercise 2: Parsing Headers and Sequences
===============================================================

Learning Objectives:
- Learn to parse FASTA format into structured data
- Separate headers from sequences
- Handle multi-line sequences
- Create data structures to store sequence information

Background:
In Exercise 1, we learned to read FASTA files and identify headers. Now we'll learn
to properly parse the format, separating each sequence from its header and handling
sequences that span multiple lines.

A typical FASTA entry looks like:
>sequence_id description text here
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC
CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC
CTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGAC

Instructions:
Complete the functions below to parse FASTA files into structured data.
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
    # TODO: Parse the FASTA content
    # Hint: Split into lines, identify headers vs sequence lines,
    #       and build sequences by concatenating sequence lines
    sequences = []
    # Your code here
    pass


def parse_fasta_dict(fasta_content):
    """
    Parse FASTA content into a dictionary with sequence IDs as keys.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        dict: Dictionary where keys are sequence IDs (without '>') 
              and values are dictionaries with 'description' and 'sequence' keys
    """
    # TODO: Parse into a more structured format
    # Example result:
    # {
    #   'sequence_1': {
    #     'description': 'Homo sapiens insulin gene',
    #     'sequence': 'ATGGCCCTGTGGATGCGC...'
    #   }
    # }
    result = {}
    # Your code here
    pass


def get_sequence_by_id(sequences_dict, sequence_id):
    """
    Retrieve a specific sequence by its ID.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        sequence_id (str): ID of the sequence to retrieve
        
    Returns:
        dict: Dictionary with sequence information, or None if not found
    """
    # TODO: Return the sequence information for the given ID
    pass


def get_all_sequence_ids(sequences_dict):
    """
    Get a list of all sequence IDs in the dataset.
    
    Args:
        sequences_dict (dict): Dictionary from parse_fasta_dict()
        
    Returns:
        list: List of all sequence IDs
    """
    # TODO: Return a list of all sequence IDs (keys)
    pass


def main():
    """Main function to test your implementations."""
    print("FASTA File Analysis - Exercise 2")
    print("=" * 40)
    
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