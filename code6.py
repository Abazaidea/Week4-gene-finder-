import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_rbs(sequence, start, upstream_range, rbs_sequence):
    upstream_seq = sequence[max(0, start - upstream_range):start]
    return rbs_sequence in upstream_seq

def find_orfs(sequence, min_length, upstream_range, rbs_sequence):
    orfs = []
    seq_len = len(sequence)
    
    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for start in range(frame, seq_len, 3):
                if seq[start:start+3] == 'ATG':
                    for end in range(start + 3, seq_len, 3):
                        if seq[end:end+3] in ['TAA', 'TAG', 'TGA']:
                            orf = seq[start:end+3]
                            if len(orf) / 3 >= min_length:
                                if strand == 1:
                                    has_rbs = find_rbs(sequence, start, upstream_range, rbs_sequence)
                                else:
                                    has_rbs = find_rbs(sequence.reverse_complement(), seq_len - end - 3, upstream_range, rbs_sequence)
                                if has_rbs:
                                    orfs.append((start, end+3, strand, orf))
                            break
    return orfs

def main(input_file, min_length, upstream_range, rbs_sequence):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = record.seq
        orfs = find_orfs(sequence, min_length, upstream_range, rbs_sequence)
        
        for start, end, strand, orf in orfs:
            strand_symbol = '+' if strand == 1 else '-'
            print(f">ORF_{start}_{end}_{strand_symbol}")
            print(orf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file with length and RBS filtering")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("-l", "--min_length", type=int, default=100, 
                        help="Minimum ORF length in codons (default: 100)")
    parser.add_argument("-u", "--upstream_range", type=int, default=20,
                        help="Number of base pairs upstream to search for RBS (default: 20)")
    parser.add_argument("-r", "--rbs_sequence", type=str, default="AGGAGG",
                        help="Ribosome Binding Site sequence to search for (default: AGGAGG)")
    args = parser.parse_args()
    
    main(args.input_file, args.min_length, args.upstream_range, args.rbs_sequence)
