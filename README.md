# Bioinformatics Algorithms & Computational Biology

This repository contains a collection of Python scripts and algorithms solving classic and advanced computational biology problems. It demonstrates proficiency in DNA/RNA sequence analysis, motif finding algorithms, and REST API integration with global biological databases like UniProt.

## 🧬 Repository Contents

### 1. Advanced Motif Finding
* **`branch_and_bound_median_search.py`**: A complete, from-scratch implementation of the **Branch and Bound Median Search** algorithm. It navigates a lexicographical tree using `NextVertex` and `ByPass` techniques to find the optimal $l$-mer motif across multiple DNA sequences by minimizing the total Hamming distance.

### 2. Classic Bioinformatics Challenges (Rosalind)
* **`rosalind_bioinformatics_challenges.py`**: An Object-Oriented toolkit (`BioinformaticsToolkit`) solving fundamental computational biology problems:
  * **Reverse Complement**: Computes the reverse complement of DNA strings.
  * **Fibonacci Rabbits**: Solves the population dynamics problem using dynamic programming.
  * **FASTA GC Content**: Parses FASTA formatted data to calculate and identify the sequence with the highest GC content.
  * **Protein Motif Search (UniProt API)**: Fetches raw protein sequences directly from the **UniProt REST API** and uses Regular Expressions (Regex) to locate overlapping N-glycosylation motifs (`N{P}[ST]{P}`).

### 3. Sequence Fundamentals
* **`dna_sequence_basics.py`**: Foundational OOP scripts for biological string validation, nucleotide counting, and DNA-to-RNA transcription.

## 🛠️ Skills & Technologies Highlighted
* **Algorithmic Complexity**: Branch and Bound, Tree Traversal, Dynamic Programming.
* **Data Parsing**: FASTA file parsing and string manipulation.
* **Web Integration**: Fetching and decoding data from the `rest.uniprot.org` API using `urllib`.
* **Pattern Matching**: Advanced Regular Expressions (`re`) with lookahead assertions for overlapping biological motifs.
