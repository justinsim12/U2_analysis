from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from typing import List


class SeqAnalyzer:

    def __init__(self, filepath: str, filetype: str):
        self.parsed_sequence = SeqIO.parse(filepath, filetype)

    def print_all_sequences(self) -> None:
        for seq_record in self.parsed_sequence:
            print(seq_record.id)
            print(repr(seq_record.seq))
            print(len(seq_record))

    def get_sequence_by_index(self, index: int) -> Seq:
        try:
            found_seq = None
            for i, seq in enumerate(self.parsed_sequence):
                if i == index:
                    print(f"Retrieved sequence: {seq.id}, at index: {index}")
                    found_seq = seq
                    break
        except Exception as e:
            print(f"Unable to retrieve sequence at index: {index}, error: {e}")
        finally:
            return found_seq

    def get_sequence_by_id(self, id: str) -> Seq:
        try:
            found_seq = None
            for i, seq in enumerate(self.parsed_sequence):
                if seq.id == id:
                    print(f"Retrieved sequence: {seq.id}, at index: {i}")
                    found_seq = seq
                    break
        except Exception as e:
            print(f"Unable to retrieve sequence with id: {id}, error: {e}")
        finally:
            return found_seq

    # very computationally complex, On^2
    def get_sequences_by_subsequence(self, sub: str | Seq | MutableSeq) -> List[Seq]:
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

