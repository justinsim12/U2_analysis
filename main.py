from SeqAnalyzer import *

analyzer = SeqAnalyzer("Jalview_U2_alignment.fa", "fasta")



if __name__ == '__main__':
    #analyzer.print_all_sequences()
    analyzer.get_sequence_by_index(0)
    analyzer.get_sequence_by_id('hs1_t2tRepeatMasker_U2#snRNA/51-120')
    analyzer.get_sequences_by_subsequence('ATTGCTTCTCGGCCTTTTGGCTAAGATCAAGCGTAAACCACAAACAGA')
