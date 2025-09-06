"""
FASTA File Analysis - Exercise 1: Basic File Reading and FASTA Format
=====================================================================

Learning Objectives:
- Understand the FASTA file format structure
- Learn how to open and read files in Python
- Practice basic string operations

Background:
FASTA format is a text-based format for representing nucleotide or protein sequences.
Each sequence has two parts:
1. Header line: starts with '>' followed by sequence identifier and description
2. Sequence lines: the actual sequence data (can span multiple lines)

Example:
>sequence_1 Human insulin gene
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC
CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC

Instructions:
Complete the functions below to analyze FASTA files.
"""

def read_fasta_file(filename):
    """
    Read a FASTA file and return its contents as a string.
    
    Args:
        filename (str): Path to the FASTA file
        
    Returns:
        str: Contents of the file
    """
    # TODO: Open and read the file, return its contents
    # Hint: Use open() function with 'r' mode
    pass


def count_sequences(fasta_content):
    """
    Count the number of sequences in FASTA content by counting header lines.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        int: Number of sequences (header lines starting with '>')
    """
    # TODO: Count lines that start with '>'
    # Hint: Split content into lines and check each line
    pass


def get_sequence_headers(fasta_content):
    """
    Extract all sequence headers from FASTA content.
    
    Args:
        fasta_content (str): Content of a FASTA file
        
    Returns:
        list: List of header lines (including the '>' character)
    """
    # TODO: Find and return all lines that start with '>'
    pass


def main():
    """Main function to test your implementations."""
    print("FASTA File Analysis - Exercise 1")
    print("=" * 40)
    
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