# U2_analysis

The goal of this project will be to hierarchically cluster the sequences and bin them into based on longest run of left to right matches.

Eg. given the following strings
aaaaaaaaabbbbbb
aaaadsvvvvvv
aaaaklkdddd
aaaaaaamsldkmflk

cluster 1:
aaaadsvvvvvv
aaaaklkdddd

cluster 2:
aaaaaaamsldkmflk

cluster 3:
aaaaaaaaabbbbbb

The data is in FASTA format, use biopython to extract the sequence records and probably best to keep the data as a seq record to store both the header and sequence.
