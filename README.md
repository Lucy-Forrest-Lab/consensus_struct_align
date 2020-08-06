# Consensus.pl 		  
Program to fuse alignments and 
calculate confidence values	  
### date: 05.05.2015 by Marcus Stamm in the Forrest Lab		  

This is a script that can fuse 4 alignments and calculate confidence values for each sequence position. 
For this purpose, 4 pairwise sequence alignments in fasta format have to be submitted.
The order of the sequences in the fasta files has to be the same in each file (i.e. reference then target). 
The first sequence is taken as the reference sequence. 
Gaps that were present in the original alignments are represented by "-".
Gaps that had to be inserted during the creation of a fused alignment are represented by a ".". 
The higher the value of the confidence score, the more reliable that column of the alignment.

### How to install
- extract the `Consensus.pl` and `FuseAlignments.pm` into the same folder

### How to use
`perl Consensus.pl -i <file1> <file2> <file3> <file4> -m <id1> <id2> <id3> <id4>  -o <outputfile>`
- the flag "-i" has to be followed by 4 filenames of alignments in fasta format 
- the flat "-m" is optional and places identifiers for the inputs into the output alignment; so they should be provided in the same order as the input files.
- the flag "-o" is optional - if no output filename is given, the results will be reported to standard output.

### Example 
- extract the folder "examples" as a folder to the same directory to which you have also extracted the script
- enter the following command: 
`perl Consensus.pl -i  ./examples/dalilite.fa ./examples/frtmalign.fa ./examples/fatcat.fa ./examples/matt.fa -m DALI FR-TM-align FATCAT MATT -o output_consensus.txt`

--

Help with converting structural alignments into fasta format: 
Please note that there are scripts to convert structure alignment outputs from DALILITE, FR-TMalign, FATCAT and MATT into fasta format in the directory /struct2fasta. 
