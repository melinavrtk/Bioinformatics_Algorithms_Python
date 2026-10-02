# =========================================================
# Chapter: Basic Bioinformatics - DNA Operations
# Features: Nucleotide Counting & DNA to RNA Transcription
# =========================================================

class DNASequence:
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        if not self.is_valid():
            raise ValueError("Invalid DNA string. Only 'A', 'C', 'G', 'T' are allowed.")

    def is_valid(self):
        """ Validates that the sequence contains only standard nucleotides. """
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        return all(nucleotide in valid_nucleotides for nucleotide in self.sequence)

    def count_nucleotides(self):
        """ Returns the count of each nucleotide in the format: A C G T """
        count_A = self.sequence.count('A')
        count_C = self.sequence.count('C')
        count_G = self.sequence.count('G')
        count_T = self.sequence.count('T')
        return count_A, count_C, count_G, count_T

    def transcribe_to_rna(self):
        """ Transcribes DNA sequence to RNA by replacing Thymine (T) with Uracil (U). """
        return self.sequence.replace('T', 'U')

# =========================================================
# MAIN EXECUTION
# =========================================================
if __name__ == "__main__":
    valid_input = False
    
    # Input Loop mimicking MATLAB logic
    while not valid_input:
        user_input = input("Sample Dataset (DNA String): ").strip().upper()
        try:
            dna = DNASequence(user_input)
            valid_input = True
        except ValueError as e:
            print(f"Error: {e}")

    # 1. Nucleotide Counting
    counts = dna.count_nucleotides()
    print(f"Nucleotide Counts (A C G T): {counts[0]} {counts[1]} {counts[2]} {counts[3]}")

    # 2. DNA to RNA Transcription
    rna_sequence = dna.transcribe_to_rna()
    print(f"Transcribed RNA: {rna_sequence}")
