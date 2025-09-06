"""
FASTA File Analysis - Exercise 4 Solution: Advanced Sequence Analysis
=====================================================================
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
    positions = []
    sequence = sequence.upper()
    pattern = pattern.upper()
    
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            positions.append(i)
    
    return positions


def translate_dna(dna_sequence, frame=0):
    """
    Translate DNA sequence to protein sequence using the genetic code.
    
    Args:
        dna_sequence (str): DNA sequence to translate
        frame (int): Reading frame (0, 1, or 2)
        
    Returns:
        str: Protein sequence (single letter amino acid codes)
    """
    dna_sequence = dna_sequence.upper()
    protein = ""
    
    # Start from the specified frame
    for i in range(frame, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3 and codon in GENETIC_CODE:
            protein += GENETIC_CODE[codon]
        else:
            protein += 'X'  # Unknown amino acid
    
    return protein


def find_start_codons(sequence):
    """
    Find all positions of start codons (ATG) in a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        list: List of positions of start codons
    """
    return find_pattern(sequence, 'ATG')


def find_stop_codons(sequence):
    """
    Find all positions of stop codons (TAA, TAG, TGA) in a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        list: List of tuples (position, codon) for stop codons
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    positions = []
    
    for stop_codon in stop_codons:
        for pos in find_pattern(sequence, stop_codon):
            positions.append((pos, stop_codon))
    
    # Sort by position
    positions.sort(key=lambda x: x[0])
    
    return positions


def find_orfs(sequence, min_length=30):
    """
    Find Open Reading Frames (ORFs) in a DNA sequence.
    
    Args:
        sequence (str): DNA sequence
        min_length (int): Minimum ORF length in nucleotides
        
    Returns:
        list: List of dictionaries with ORF information
    """
    orfs = []
    sequence = sequence.upper()
    
    # Find all start positions
    start_positions = find_start_codons(sequence)
    
    for start_pos in start_positions:
        # Look for stop codon in the same reading frame
        for i in range(start_pos + 3, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in ['TAA', 'TAG', 'TGA']:
                # Found stop codon
                orf_length = i + 3 - start_pos
                if orf_length >= min_length:
                    orf_sequence = sequence[start_pos:i+3]
                    protein = translate_dna(orf_sequence)
                    
                    orfs.append({
                        'start': start_pos,
                        'end': i + 3,
                        'length': orf_length,
                        'protein': protein
                    })
                break
    
    return orfs


def reverse_complement(sequence):
    """
    Calculate the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: Reverse complement sequence
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    sequence = sequence.upper()
    
    # Create complement
    complement = ""
    for nucleotide in sequence:
        if nucleotide in complement_map:
            complement += complement_map[nucleotide]
        else:
            complement += nucleotide  # Keep unknown nucleotides
    
    # Reverse the complement
    return complement[::-1]


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
        '*': 0,  # Stop codon
        'X': 110  # Average amino acid weight for unknown
    }
    
    total_weight = 0
    protein_sequence = protein_sequence.upper()
    
    for amino_acid in protein_sequence:
        if amino_acid in aa_weights:
            total_weight += aa_weights[amino_acid]
        else:
            total_weight += aa_weights['X']  # Use average for unknown
    
    return total_weight


def main():
    """Main function to test the implementations."""
    print("FASTA File Analysis - Exercise 4 Solution")
    print("=" * 45)
    
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
            print(f"Stop codons: {stop_codons[:5]}")  # Show first 5
            
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
            print(f"Original (first 50 nt): {sequence[:50]}")
            print(f"Reverse complement: {rev_comp}")
            
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