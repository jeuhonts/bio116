"""
FASTA File Analysis - Exercise 1 Solution: Basic File Reading and FASTA Format
==============================================================================
"""

def read_fasta_file(filename):
    """
    Read a FASTA file and return its contents as a string.
    
    Args:
        filename (str): Path to the FASTA file
        
    Returns:
        str: Contents of the file
    """
    with open(filename, 'r') as file:
        return file.read()


def count_sequences(fasta_content):
    """
    Count the number of sequences in FASTA content by counting header lines.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        int: Number of sequences (header lines starting with '>')
    """
    lines = fasta_content.split('\n')
    count = 0
    for line in lines:
        if line.startswith('>'):
            count += 1
    return count


def get_sequence_headers(fasta_content):
    """
    Extract all sequence headers from FASTA content.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        list: List of header lines (including the '>' character)
    """
    lines = fasta_content.split('\n')
    headers = []
    for line in lines:
        if line.startswith('>'):
            headers.append(line)
    return headers


def main():
    """Main function to test the implementations."""
    print("FASTA File Analysis - Exercise 1 Solution")
    print("=" * 45)
    
    # Test with sample sequences
    filename = "data/sample_sequences.fasta"
    
    try:
        # Read the file
        content = read_fasta_file(filename)
        print(f"Successfully read file: {filename}")
        print(f"File size: {len(content)} characters")
        print()
        
        # Count sequences
        seq_count = count_sequences(content)
        print(f"Number of sequences found: {seq_count}")
        print()
        
        # Get headers
        headers = get_sequence_headers(content)
        print("Sequence headers:")
        for i, header in enumerate(headers, 1):
            print(f"  {i}. {header}")
        print()
        
        # Show first 100 characters of the file
        print("First 100 characters of the file:")
        print(content[:100] + "..." if len(content) > 100 else content)
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'")
        print("Make sure you're running this script from the correct directory.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()