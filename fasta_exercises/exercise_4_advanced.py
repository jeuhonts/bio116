"""
FASTA File Analysis - Exercise 4: Advanced Sequence Analysis
============================================================

Learning Objectives:
- Find patterns and motifs in sequences
- Translate DNA to protein sequences
- Search for open reading frames (ORFs)
- Calculate molecular properties
- Work with the genetic code

Background:
Advanced sequence analysis includes finding specific patterns, translating DNA to proteins,
and identifying potential coding regions. The genetic code is used to translate triplets
of nucleotides (codons) into amino acids.

Key concepts:
- Codon: A sequence of three nucleotides that codes for an amino acid
- Start codon: ATG (codes for Methionine, starts protein synthesis)
- Stop codons: TAA, TAG, TGA (terminate protein synthesis)
- Open Reading Frame (ORF): A sequence that starts with a start codon and ends with a stop codon

Instructions:
Complete the functions below for advanced sequence analysis.
"""

# Standard genetic code dictionary
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def find_pattern(sequence, pattern):
    """
    Find all occurrences of a pattern in a sequence.
    
    Args:
        sequence (str): DNA sequence to search
        pattern (str): Pattern to find
        
    Returns:
        list: List of positions where pattern starts (0-based indexing)
    """
    # TODO: Find all starting positions of the pattern
    positions = []
    # Your code here
    pass


def translate_dna(dna_sequence, frame=0):
    """
    Translate DNA sequence to protein sequence using the genetic code.
    
    Args:
        dna_sequence (str): DNA sequence to translate
        frame (int): Reading frame (0, 1, or 2)
        
    Returns:
        str: Protein sequence (single letter amino acid codes)
    """
    # TODO: Translate DNA to protein
    # Start from the specified frame
    # Group nucleotides into codons (triplets)
    # Look up each codon in the genetic code
    protein = ""
    # Your code here
    pass


def find_start_codons(sequence):
    """
    Find all positions of start codons (ATG) in a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        list: List of positions of start codons
    """
    # TODO: Find all ATG positions
    pass


def find_stop_codons(sequence):
    """
    Find all positions of stop codons (TAA, TAG, TGA) in a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        list: List of tuples (position, codon) for stop codons
    """
    # TODO: Find all stop codon positions
    stop_codons = ['TAA', 'TAG', 'TGA']
    positions = []
    # Your code here
    pass


def find_orfs(sequence, min_length=30):
    """
    Find Open Reading Frames (ORFs) in a DNA sequence.
    An ORF starts with ATG and ends with a stop codon.
    
    Args:
        sequence (str): DNA sequence
        min_length (int): Minimum ORF length in nucleotides
        
    Returns:
        list: List of dictionaries with ORF information:
              {'start': position, 'end': position, 'length': length, 'protein': protein_sequence}
    """
    # TODO: Find all ORFs
    # Look for start codons, then find the next stop codon in frame
    # Calculate length and translate to protein
    orfs = []
    # Your code here
    pass


def reverse_complement(sequence):
    """
    Calculate the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: Reverse complement sequence
    """
    # TODO: Create reverse complement
    # A↔T, G↔C, reverse the order
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    # Your code here
    pass


def calculate_molecular_weight(protein_sequence):
    """
    Calculate approximate molecular weight of a protein sequence.
    
    Args:
        protein_sequence (str): Protein sequence (amino acid codes)
        
    Returns:
        float: Molecular weight in Daltons
    """
    # Approximate molecular weights of amino acids (in Daltons)
    aa_weights = {
        'A': 89, 'C': 121, 'D': 133, 'E': 147, 'F': 165,
        'G': 75, 'H': 155, 'I': 131, 'K': 146, 'L': 131,
        'M': 149, 'N': 132, 'P': 115, 'Q': 146, 'R': 174,
        'S': 105, 'T': 119, 'V': 117, 'W': 204, 'Y': 181,
        '*': 0  # Stop codon
    }
    
    # TODO: Calculate total molecular weight
    total_weight = 0
    # Your code here
    pass


def main():
    """Main function to test your implementations."""
    print("FASTA File Analysis - Exercise 4")
    print("=" * 40)
    
    filename = "data/sample_sequences.fasta"
    
    try:
        # Read and parse the file
        with open(filename, 'r') as file:
            content = file.read()
        
        # Simple parsing function
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
        
        sequences_list = parse_fasta_simple(content)
        
        print(f"Advanced analysis of sequences from: {filename}")
        print()
        
        # Analyze first sequence in detail
        if sequences_list:
            header, sequence = sequences_list[0]
            seq_id = header.split()[0][1:]
            
            print(f"Detailed analysis of: {seq_id}")
            print(f"Sequence length: {len(sequence)}")
            print()
            
            # Find patterns
            atg_positions = find_pattern(sequence, 'ATG')
            print(f"ATG positions: {atg_positions}")
            
            # Find start and stop codons
            start_codons = find_start_codons(sequence)
            stop_codons = find_stop_codons(sequence)
            print(f"Start codons (ATG): {start_codons}")
            print(f"Stop codons: {stop_codons}")
            
            # Translation in different frames
            for frame in range(3):
                protein = translate_dna(sequence, frame)
                print(f"Frame {frame} translation (first 20 AA): {protein[:20]}")
            
            # Find ORFs
            orfs = find_orfs(sequence)
            print(f"Found {len(orfs)} ORFs (min 30 nt)")
            for i, orf in enumerate(orfs[:3]):  # Show first 3 ORFs
                print(f"  ORF {i+1}: pos {orf['start']}-{orf['end']}, "
                      f"length {orf['length']}, protein: {orf['protein'][:10]}...")
            
            # Reverse complement
            rev_comp = reverse_complement(sequence[:50])  # First 50 nt
            print(f"Reverse complement (first 50 nt): {rev_comp}")
            
            # Molecular weight calculation
            if orfs:
                first_orf_protein = orfs[0]['protein']
                mol_weight = calculate_molecular_weight(first_orf_protein)
                print(f"Molecular weight of first ORF: {mol_weight:.2f} Da")
            
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'")
        print("Make sure you're running this script from the correct directory.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()