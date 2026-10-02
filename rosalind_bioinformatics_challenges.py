# =========================================================
# Chapter: Classic Bioinformatics Algorithms (Rosalind Challenges)
# Features: Reverse Complement, Fibonacci, GC Content, UniProt API
# =========================================================
import urllib.request
import re

class BioinformaticsToolkit:
    
    @staticmethod
    def reverse_complement(dna_sequence):
        """
        Problem 1: Finds the reverse complement of a DNA string.
        (A -> T, T -> A, C -> G, G -> C)
        """
        # Dictionary for complementarity
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # Reverse the sequence and map each nucleotide
        rev_comp = "".join(complement_map[base] for base in reversed(dna_sequence))
        return rev_comp

    @staticmethod
    def fibonacci_rabbits(months, offspring_pairs):
        """
        Problem 2: Wascally Wabbits (Fibonacci Sequence adaptation).
        Calculates total rabbit pairs after 'n' months, where each 
        mature pair produces 'k' offspring pairs.
        """
        if months == 1 or months == 2:
            return 1
            
        F = [0] * (months + 1)
        F[1], F[2] = 1, 1
        
        for i in range(3, months + 1):
            F[i] = F[i-1] + offspring_pairs * F[i-2]
            
        return F[months]

    @staticmethod
    def highest_gc_content(fasta_content):
        """
        Problem 3: Computes GC Content from a FASTA formatted string
        and returns the ID with the highest GC percentage.
        """
        lines = fasta_content.strip().split('\n')
        max_gc = 0.0
        max_id = ""
        
        current_id = ""
        current_seq = ""
        
        def calculate_gc(sequence):
            if not sequence: return 0.0
            gc_count = sequence.count('G') + sequence.count('C')
            return (gc_count / len(sequence)) * 100

        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Process the previous sequence before starting the new one
                if current_seq:
                    gc = calculate_gc(current_seq)
                    if gc > max_gc:
                        max_gc = gc
                        max_id = current_id
                        
                current_id = line[1:] # Remove '>'
                current_seq = ""
            else:
                current_seq += line

        # Process the final sequence in the file
        if current_seq:
            gc = calculate_gc(current_seq)
            if gc > max_gc:
                max_gc = gc
                max_id = current_id
                
        return max_id, max_gc

    @staticmethod
    def find_protein_motif(uniprot_ids):
        """
        Problem 4 (Problem 06 from Assignment): Finding a Protein Motif.
        Fetches FASTA sequences from UniProt API and finds N-glycosylation motifs.
        Motif Regex: N{P}[ST]{P} -> N followed by not P, followed by S or T, followed by not P.
        """
        results = {}
        # Regex for N-glycosylation: N, followed by anything but P, followed by S or T, followed by anything but P
        motif_regex = re.compile(r'(?=(N[^P][ST][^P]))')
        
        for protein_id in uniprot_ids:
            protein_id = protein_id.strip()
            url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
            
            try:
                # Fetch data directly from UniProt API
                response = urllib.request.urlopen(url)
                fasta_data = response.read().decode('utf-8')
                
                # Split and extract only the sequence (ignore header line)
                lines = fasta_data.strip().split('\n')
                sequence = "".join(line.strip() for line in lines[1:])
                
                # Find overlapping matches using lookahead regex
                # Note: Biological indexing starts at 1, not 0
                positions = [match.start() + 1 for match in motif_regex.finditer(sequence)]
                
                if positions:
                    results[protein_id] = positions
                    
            except Exception as e:
                print(f"Failed to fetch or process UniProt ID {protein_id}: {e}")
                
        return results

# =========================================================
# MAIN EXECUTION & TESTING
# =========================================================
if __name__ == "__main__":
    toolkit = BioinformaticsToolkit()
    
    print("--- 1. Reverse Complement ---")
    dna_sample = "GTCA"
    print(f"Original: {dna_sample}")
    print(f"Reverse Complement: {toolkit.reverse_complement(dna_sample)}")
    
    print("\n--- 2. Fibonacci Rabbits ---")
    months (n), offspring (k) = 5, 3
    print(f"Total pairs after {n} months (k={k}): {toolkit.fibonacci_rabbits(n, k)}")
    
    print("\n--- 3. Highest GC Content ---")
    sample_fasta = ">Rosalind_6404\nCCTGCGGAAGATCGGCACTCGGTGTCA\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTT\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCT"
    best_id, best_gc = toolkit.highest_gc_content(sample_fasta)
    print(f"ID: {best_id}\nGC Content: {best_gc:.6f}%")
    
    print("\n--- 4. Protein Motif Search (UniProt API) ---")
    test_ids = ["B5ZC00", "P07204", "P20840"]
    motif_results = toolkit.find_protein_motif(test_ids)
    
    for pid, pos_list in motif_results.items():
        pos_str = " ".join(map(str, pos_list))
        print(f"{pid}\n{pos_str}")
