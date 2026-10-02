# =========================================================
# Chapter: Motif Finding Algorithms
# Algorithm: Branch and Bound Median Search
# =========================================================
import numpy as np

def letter_to_number(dna_strings):
    """ Converts DNA strings ('A', 'C', 'G', 'T') to numeric representations (1, 2, 3, 4). """
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return [[mapping[char] for char in seq] for seq in dna_strings]

def number_to_letter(numeric_array):
    """ Converts numeric motif representation back to a DNA string. """
    mapping = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
    return "".join([mapping[num] for num in numeric_array])

def next_vertex(a, i, L, k):
    """ Generates the next node in the lexicographical tree. """
    if i < k:
        a[i] = 1
        return a, i + 1
    else:
        for j in range(k - 1, -1, -1):
            if a[j] < L:
                a[j] += 1
                return a, j + 1
    return a, 0

def bypass(a, i, L, k):
    """ Bypasses the subtree if the current prefix is unpromising. """
    for j in range(i - 1, -1, -1):
        if a[j] < L:
            a[j] += 1
            return a, j + 1
    return a, 0

def prefix_total_distance(prefix, dna_matrix, prefix_length):
    """ 
    Computes the total Hamming distance of a motif prefix against all DNA sequences.
    This provides the optimistic lower bound for the B&B search.
    """
    t = len(dna_matrix)
    n = len(dna_matrix[0])
    total_dist = 0
    
    for row in range(t):
        min_dist = float('inf')
        for col in range(n - prefix_length + 1):
            # Calculate Hamming distance between prefix and DNA substring
            dist = sum(1 for x, y in zip(prefix, dna_matrix[row][col : col + prefix_length]) if x != y)
            if dist < min_dist:
                min_dist = dist
        total_dist += min_dist
        
    return total_dist

def branch_and_bound_median_search(dna_matrix, t, n, k):
    """
    Implements the Branch and Bound Median Search algorithm for motif finding.
    """
    L = 4  # Alphabet size (A, C, G, T)
    dna_nums = letter_to_number(dna_matrix)
    
    # Initialization
    s = [1] * k
    best_distance = float('inf')
    best_word = s.copy()
    
    i = 1
    while i > 0:
        if i < k:
            # Check the optimistic distance of the prefix (Bounding)
            optimistic_dist = prefix_total_distance(s[:i], dna_nums, i)
            if optimistic_dist > best_distance:
                s, i = bypass(s, i, L, k)  # Prune the tree
            else:
                s, i = next_vertex(s, i, L, k) # Dig deeper
        else:
            # Leaf node reached: evaluate full motif
            current_dist = prefix_total_distance(s, dna_nums, k)
            if current_dist < best_distance:
                best_distance = current_dist
                best_word = s.copy()
            s, i = next_vertex(s, i, L, k)
            
    return number_to_letter(best_word), best_distance

# =========================================================
# MAIN EXECUTION (Running Example)
# =========================================================
if __name__ == "__main__":
    # Test Data based on assignment specifications
    DNA = [
        'ATCCAGCT',
        'GGGCAACT',
        'ATGGATCT',
        'AAGCAACC',
        'TTGGAACT',
        'ATGCCATT',
        'ATGGCACT'
    ]
    
    t = len(DNA)
    n = len(DNA[0])
    lmer = 7 # Motif length
    
    print("Executing Branch and Bound Median Search...")
    print(f"Finding {lmer}-mer motif in {t} sequences of length {n}...\n")
    
    best_motif, best_dist = branch_and_bound_median_search(DNA, t, n, lmer)
    
    print("--- RESULTS ---")
    print(f"Best Distance (BD): {best_dist}")
    print(f"Best Motif (BW): {best_motif}")
