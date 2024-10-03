# BioE 230 Week 3-4

a simple tool that finds genes in genome sequence using Python
- I used chatgpt to help with this HW.
- I also got help form my classmates Ashhad and Omar.
- my prompt to chatgpt beacuse I had not functioning code, I wote the code and asked chat gpt to fix the code to get the final version to be listed below
#### Chatgpt Promt
we will write a simple tool that finds genes in genome sequence using Python. The tool you write needs to take command line parameters
for the input file. The input file should consist of a single FASTA file
containing a single genome.
Your tool should read a FASTA file, and output any region between a
start (‘ATG’) and stop codon (‘TAA’, ‘TAG’, ‘TGA’) in that FASTA file;
you must consider to include the reverse complement and search
in all six possible reading frames for genes. use BioPython. 
The of the code output should be in this form 

```javascript
 >ORF_#
AGGAGATCCGATCGATAGCTGAATAGATAGCTGAGATAGCAGCTGCTCTCTTTTCCCAAA"
```

## Question 1 

#### The tool you write needs to take command line parameters for the input file. The input file should consist of a single FASTA file containing a single genome. Your tool should read a FASTA file, and output any region between a start (‘ATG’) and stop codon (‘TAA’, ‘TAG’, ‘TGA’) in that FASTA file; you must consider three possible reading frames but may ignore reverse compliments. 

### Command 
```javascript
python code1.py /home/abazaiea/ncbi_dataset/data/GCF_000008725.1_ASM872v1_genomic.fna > output1.txt
```


###  Output
```javascript
File output1.txt 

```
## Question 2
#### Extend your tool to include the reverse complement and search in all six possible reading frames for genes. You could reuse the code you wrote in the first week’s assignments for generating reverse compliments of a string, or use BioPython to complete the assignment.

### Command 
```javascript
python code2.py /home/abazaiea/ncbi_dataset/data/GCF_000008725.1_ASM872v1_genomic.fna > output2.txt
```


###  Output
```javascript
File output2.txt 

```

##  Question 3

#### Use your code to solve the Open Reading Frame problem on Rosalind (Problem 72)


### Command 
```javascript
python code3.py /home/abazaiea/ncbi_dataset/data/rosalind.fna > output3.txt
```


###  Output
```javascript
File output3.txt 

```
## Question 4
#### Apply your code to the genomes you downloaded and find all Open Reading Frames in the 14 genomes you downloaded. Do not run your command 14 times, write a single command to find all Open Reading Frames in all 14 genomes! (Hint: use the for command on bash!)

### command 
```javascript
find /home/abazaiea/ncbi_dataset/data -type f -name "*GCF*.fna" | while read genome; do python gene2.py "$genome"; done > output4.txt
```

### output
```javascript
File output4.txt

```

## Question 5
#### Now you have Open Reading Frames, but not all of them will be genes. Implement a filter by length: discard short ORFs that are unlikely to be functional genes (e.g., less than 100 codons, but make the length a parameter of your tool).



### command 
```javascript
python code5.py /home/abazaiea/ncbi_dataset/data/GCF_000008725.1_ASM872v1_genomic.fna -l 100 > output5.txt
```

### output
```javascript
File output5.txt

```
## Question 6

#### Even when filtered by length, not all ORFs will eventually be translated into a protein. One thing to look for is a ribosome binding site which is usually located 4-20bp upstream of the start coding. Scan upstream of the predicted start codon (e.g., -20 bp, but make this a parameter of your tool) for a ribosome binding site (RBS); there are several sequences that can indicate an RBS, the most common being the Shine-Dalgarno sequence (AGGAGG). Filter all predicted ORFs based on whether they contain a Shine-Dalgarno sequence up to 20bp upstream of the start codon.



### command 
```javascript
python code6.py /home/abazaiea/ncbi_dataset/data/GCF_000008725.1_ASM872v1_genomic.fna -l 100 > output6.txt

```

### output
```javascript
File output6.txt

```
## Final Gene Finder Tool 

```javascript
import argparse  # Importing argparse to handle command-line arguments
from Bio import SeqIO  # Importing SeqIO from Biopython to read sequence files
from Bio.Seq import Seq  # Importing Seq class from Biopython for sequence operations

def find_rbs(sequence, start, upstream_range, rbs_sequence):
    """
    Check if the specified RBS sequence is present upstream of a given start position in the sequence.
    
    Args:
        sequence: The nucleotide sequence to search within.
        start: The position of the start codon.
        upstream_range: The number of bases to check upstream of the start codon.
        rbs_sequence: The RBS sequence to look for.
    
    Returns:
        True if the RBS sequence is found; False otherwise.
    """
    # Get the sequence segment upstream of the start codon
    upstream_seq = sequence[max(0, start - upstream_range):start]
    # Check if the RBS sequence is in the upstream sequence
    return rbs_sequence in upstream_seq

def find_orfs(sequence, min_length, upstream_range, rbs_sequence):
    """
    Find Open Reading Frames (ORFs) in the given sequence.
    
    Args:
        sequence: The nucleotide sequence to search for ORFs.
        min_length: Minimum length of the ORFs (in codons).
        upstream_range: The range to check for the RBS upstream of the start codon.
        rbs_sequence: The RBS sequence to search for.
    
    Returns:
        A list of found ORFs with their start, end positions, strand information, and sequence.
    """
    orfs = []  # Initialize an empty list to store ORFs
    seq_len = len(sequence)  # Get the length of the sequence
    
    # Check both strands (forward and reverse complement)
    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        # Check all three reading frames
        for frame in range(3):
            # Iterate over possible start positions for codons
            for start in range(frame, seq_len, 3):
                # Check if the start codon (ATG) is found
                if seq[start:start+3] == 'ATG':
                    # Iterate through possible stop codons
                    for end in range(start + 3, seq_len, 3):
                        if seq[end:end+3] in ['TAA', 'TAG', 'TGA']:
                            # Get the ORF sequence from start to end
                            orf = seq[start:end+3]
                            # Check if the ORF length is greater than or equal to min_length
                            if len(orf) / 3 >= min_length:
                                # Check for RBS presence based on the strand
                                if strand == 1:
                                    has_rbs = find_rbs(sequence, start, upstream_range, rbs_sequence)
                                else:
                                    has_rbs = find_rbs(sequence.reverse_complement(), seq_len - end - 3, upstream_range, rbs_sequence)
                                
                                # If RBS is found, add the ORF to the list
                                if has_rbs:
                                    orfs.append((start, end+3, strand, orf))
                            break  # Stop looking for further stop codons after finding one
    return orfs  # Return the list of found ORFs

def main(input_file, min_length, upstream_range, rbs_sequence):
    """
    Main function to process the input FASTA file and find ORFs.
    
    Args:
        input_file: Path to the input FASTA file.
        min_length: Minimum ORF length in codons.
        upstream_range: Number of base pairs upstream to search for RBS.
        rbs_sequence: Ribosome Binding Site sequence to search for.
    """
    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = record.seq  # Extract the sequence from the record
        # Find ORFs in the sequence
        orfs = find_orfs(sequence, min_length, upstream_range, rbs_sequence)
        
        # Print each found ORF with its details
        for start, end, strand, orf in orfs:
            strand_symbol = '+' if strand == 1 else '-'  # Set strand symbol based on the strand
            print(f">ORF_{start}_{end}_{strand_symbol}")  # Print header with ORF details
            print(orf)  # Print the ORF sequence

if __name__ == "__main__":
    # Set up argument parser for command-line execution
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file with length and RBS filtering")
    parser.add_argument("input_file", help="Input FASTA file")  # Required input file argument
    parser.add_argument("-l", "--min_length", type=int, default=100, 
                        help="Minimum ORF length in codons (default: 100)")  # Optional argument for min ORF length
    parser.add_argument("-u", "--upstream_range", type=int, default=20,
                        help="Number of base pairs upstream to search for RBS (default: 20)")  # Optional argument for upstream range
    parser.add_argument("-r", "--rbs_sequence", type=str, default="AGGAGG",
                        help="Ribosome Binding Site sequence to search for (default: AGGAGG)")  # Optional argument for RBS sequence
    args = parser.parse_args()  # Parse the command-line arguments
    
    # Call the main function with the parsed arguments
    main(args.input_file, args.min_length, args.upstream_range, args.rbs_sequence)


```
