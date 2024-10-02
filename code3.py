import sys
import os

# Define stop codons
STOP_CODONS = ["TAA", "TAG", "TGA"]

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ""
    
    for nucleotide in reversed(sequence):
        rev_comp += complement.get(nucleotide, 'N')  # 'N' for unknown nucleotides
    return rev_comp

def find_genes(sequence):
    """Find genes between start (ATG) and stop codons in three reading frames."""
    sequence = sequence.upper()  # Convert the sequence to uppercase
    for frame in range(3):
        start = -1
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i + 3]
            if codon == "ATG":  # Start codon
                start = i
            elif codon in STOP_CODONS and start != -1:  # Stop codon
                gene = sequence[start:i + 3]
                gene_length = len(gene) // 3  # Number of codons
                print("Gene found from position {} to {}: {} (Length: {} codons)".format(
                    start, i + 3, gene, gene_length))
                start = -1  # Reset to find the next gene

def translate_dna_to_protein(dna_sequence):
    """Translate a DNA sequence to a protein string."""
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '',   'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '',   'TGG': 'W',
    }

    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid == '':
            break  # Stop translation at the stop codon
        protein += amino_acid
    return protein

def find_proteins(sequence):
    """Find all distinct candidate protein strings from ORFs."""
    proteins = set()
    
    # Check all six reading frames
    for frame in range(6):
        seq = sequence[frame % 3:] if frame < 3 else reverse_complement(sequence)[frame % 3:]
        start_positions = [i for i in range(len(seq)) if seq.startswith("ATG", i)]
        
        for start in start_positions:
            orf = ""
            for i in range(start, len(seq) - 2, 3):
                codon = seq[i:i + 3]
                if codon in STOP_CODONS:
                    if orf:  # Only add protein if ORF is not empty
                        protein = translate_dna_to_protein(orf)
                        if protein:  # Ensure non-empty protein
                            proteins.add(protein)
                    break
                orf += codon  # Append to ORF until a stop codon is found

    return proteins

def read_fasta(filename):
    """Read a FASTA (or .fna) file and return the sequence."""
    if not os.path.isfile(filename):  # Check if the file exists
        print("Error: File '{}' does not exist.".format(filename))
        sys.exit(1)
    
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(">"):  # Ignore header lines starting with '>'
                sequence += line.strip()
    return sequence

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    # Check if the file has a proper extension
    if not fasta_file.endswith(('.fna', '.fasta', '.fa')):
        print("Error: Please provide a valid FASTA file (.fna, .fasta, .fa)")
        sys.exit(1)
    
    # Read the sequence
    genome_sequence = read_fasta(fasta_file)
    
    # Find genes in the original sequence
    print("Finding genes in the original sequence:")
    find_genes(genome_sequence)

    # Find proteins from the original sequence
    proteins = find_proteins(genome_sequence)
    print("\nDistinct proteins found in the original sequence:")
    for protein in proteins:
        print(protein)

    # Find genes in the reverse complement
    rev_genome_sequence = reverse_complement(genome_sequence)
    print("\nFinding genes in the reverse complement sequence:")
    find_genes(rev_genome_sequence)

    # Find proteins from the reverse complement
    rev_proteins = find_proteins(rev_genome_sequence)
    print("\nDistinct proteins found in the reverse complement sequence:")
    for protein in rev_proteins:
        print(protein)

if __name__ == "__main__":
    main()
