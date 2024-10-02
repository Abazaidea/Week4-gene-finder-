# BioE 230 Week 3-4

a simple tool that finds genes in genome sequence using Python

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


```
