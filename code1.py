import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence):
    orfs = []
    seq_len = len(sequence)
    
    for frame in range(3):
        for start in range(frame, seq_len, 3):
            codon = sequence[start:start+3]
            if codon == 'ATG':
                for end in range(start + 3, seq_len, 3):
                    stop_codon = sequence[end:end+3]
                    if stop_codon in ['TAA', 'TAG', 'TGA']:
                        orf = sequence[start:end+3]
                        orfs.append((start, end+3, orf))
                        break
    return orfs

def main(input_file):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        orfs = find_orfs(sequence)
        
        for start, end, orf in orfs:
            print(f">ORF_{start}_{end}")
            print(orf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()
    
    main(args.input_file)
