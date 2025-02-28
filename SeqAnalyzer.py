from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from typing import List
import re


class SeqAnalyzer:

    def __init__(self, filepath: str, filetype: str):
        self.parsed_sequence = []
        self.index_sequence_map = {}
        self.id_sequence_map = {}
        self.chromosome_range_map = {}
        for index, rec in enumerate(SeqIO.parse(filepath, filetype)):
            match = re.search(r"range=([^\s]+)", rec.description)
            range = match.group(1) if match else None
            self.chromosome_range_map[range] = rec
            self.parsed_sequence.append(rec)
            self.index_sequence_map[index] = rec
            self.id_sequence_map[rec.id] = rec
            #self.range_sequence_map

            
    def print_all_sequences(self) -> None:
        for seq_record in self.parsed_sequence:
            print(seq_record.id)
            print(repr(seq_record.seq))
            print(len(seq_record))

    def get_sequence_by_index(self, index: int) -> SeqRecord:
        seq = self.index_sequence_map.get(index, None)

        if seq:
            print(f"Retrieved sequence: {seq.id}, at index: {index}")
            return seq
        else:      
            print(f"Unable to retrieve sequence at index: {index}")
        

    def get_sequence_by_id(self, id: str) -> SeqRecord:
        seq = self.id_sequence_map.get(id, None)

        if seq:
            print(f"Retrieved sequence: {seq.id}")
            return seq
        else:      
            print(f"Unable to retrieve sequence with id: {id}")

    def get_sequence_by_chromosome_range(self, range: str) -> SeqRecord:
        seq = self.chromosome_range_map.get(range, None)

        if seq:
            print(f"Retrieved sequence: {range}")
            return seq
        else:      
            print(f"Unable to retrieve sequence with range: {range}")


    # very computationally complex, On^2
    def get_sequences_by_subsequence(self, sub: str | SeqRecord | MutableSeq) -> List[SeqRecord]:
        try:
            found_sequences = []
            for i, seq in enumerate(self.parsed_sequence):
                if seq.seq.find(sub) != -1:
                    print(f"Retrieved sequence: {seq.id}, at index: {i}")
                    found_sequences.append(seq)
        except Exception as e:
            print(f"Unable to retrieve sequences, error: {e}")
        finally:
            return found_sequences
        
    def get_longest_common_sequence(self, input1: str | int | Seq, input2: str | int | Seq) -> Seq | str:
        seq1, seq2 = None, None
        if isinstance(input1, int):
            seq1 = self.index_sequence_map.get(input1)
        if isinstance(input2, int):
            seq2 = self.index_sequence_map.get(input2)
        if isinstance(input1, Seq):
            seq1 = input1
        if isinstance(input2, Seq):
            seq2 = input2
        if isinstance(input1, str) and 'hs1' in input1:
            seq1 = self.id_sequence_map.get(input1)
        if isinstance(input1, str):
            seq1 = input1
        if isinstance(input2, str) and 'hs1' in input2:
            seq1 = self.id_sequence_map.get(input2)
        if isinstance(input2, str):
            seq1 = input2
        #common dynamic programming logic to find longest common substring
        m, n = len(seq1), len(seq2)
        dp = [[0] * (n + 1) for _ in range(m + 1)]  # DP table initialized with 0
        max_length = 0  # Stores the length of the longest common substring
        end_index = 0   # Stores the ending index of LCS in str1

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if seq1[i - 1] == seq2[j - 1]:  # Characters match
                    dp[i][j] = dp[i - 1][j - 1] + 1
                    
                    if dp[i][j] > max_length:  # Update max_length and end_index
                        max_length = dp[i][j]
                        end_index = i

        # Extract the longest common substring using end_index and max_length
        longest_subsequence = seq1[end_index - max_length:end_index]
        print(f"Retrieved longest common subsequence: {longest_subsequence}")
        return longest_subsequence



