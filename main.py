from SeqAnalyzer import *
from pprint import pprint
import re

analyzer = SeqAnalyzer("Jalview_U2_alignment.fa", "fasta")

'''
Some random chromosome ranges
chr17:22821101-22821389
chr4:86914288-86914576
chr12:9138151-9138439
chr3:68894838-68895126
chr4:3177299-3177586
chr14:74977151-74977439
chr14:24233967-24234254
chr11:102851103-102851391

'''

if __name__ == '__main__':
    #analyzer.print_all_sequences()
    #seq1 = analyzer.get_sequence_by_index(0)
    #seq2 = analyzer.get_sequence_by_index(3)

    pprint(len(analyzer.chromosome_range_map.keys()))

    seq1 = analyzer.get_sequence_by_chromosome_range('chr4:86914288-86914576')
    seq2 = analyzer.get_sequence_by_chromosome_range('chr14:74977151-74977439')

    seq1_string = str(seq1.seq)
    seq2_string = str(seq2.seq)

    print(seq1_string)
    print(seq2_string)

    print(len(seq1_string))
    print(len(seq2_string))
    
    #analyzer.get_sequences_by_subsequence('ATTGCTTCTCGGCCTTTTGGCTAAGATCAAGCGTAAACCACAAACAGA')
    analyzer.get_longest_common_sequence(seq1.seq, seq2.seq)


   

