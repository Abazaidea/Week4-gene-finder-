import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_length):
    orfs = []
    seq_len = len(sequence)
    
    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for start in range(frame, seq_len, 3):
                if seq[start:start+3] == 'ATG':
                    for end in range(start + 3, seq_len, 3):
                        if seq[end:end+3] in ['TAA', 'TAG', 'TGA']:
                            orf = seq[start:end+3]
                            if len(orf) / 3 >= min_length:  # Check if ORF meets minimum length
                                orfs.append((start, end+3, strand, orf))
                            break
    return orfs

def main(input_file, min_length):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = record.seq
        orfs = find_orfs(sequence, min_length)
        
        for start, end, strand, orf in orfs:
            strand_symbol = '+' if strand == 1 else '-'
            print(f">ORF_{start}_{end}_{strand_symbol}")
            print(orf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("-l", "--min_length", type=int, default=100, 
                        help="Minimum ORF length in codons (default: 100)")
    args = parser.parse_args()
    
    main(args.input_file, args.min_length)
