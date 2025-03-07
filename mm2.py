#!/usr/bin/env python3

import argparse

# ----------------------------
# Helper Functions for File Parsing
# ----------------------------

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequences.
    Key: sequence identifier (without the '>').
    Value: nucleotide sequence as a string.
    """
    sequences = {}
    with open(fasta_file, 'r') as fa_fh:
        seq_id = None
        seq_lines = []
        for line in fa_fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_id:
                    sequences[seq_id] = ''.join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.upper())  # Ensure consistency
        if seq_id:
            sequences[seq_id] = ''.join(seq_lines)
    return sequences

def parse_motifs(motif_file):
    """
    Parse the motif file (one motif per line) and return a list of motifs.
    """
    motifs = []
    with open(motif_file, 'r') as motif_fh:
        for line in motif_fh:
            motif = line.strip().upper()
            if motif:
                motifs.append(motif)
    return motifs

# ----------------------------
# Motif Expansion Class
# ----------------------------

class MotifBranch:
    ambiguous_map = {
        'Y': ['T', 'C'],
        'R': ['A', 'G'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'N': ['A', 'T', 'G', 'C']
    }

    def __init__(self, sequence=''):
        self.sequence = sequence

    def expand_all(self, motif_pattern):
        """
        Recursively expand an ambiguous motif pattern into all possible sequences.
        """
        if not motif_pattern:
            return {self.sequence}  # Return a set of sequences

        char = motif_pattern[0]
        next_sequences = [self.sequence + letter for letter in MotifBranch.ambiguous_map.get(char, [char])]

        results = set()
        for seq in next_sequences:
            results.update(MotifBranch(seq).expand_all(motif_pattern[1:]))
        return results

def expand_motif(motif):
    return MotifBranch().expand_all(motif)

# ----------------------------
# FASTA Parser with Sliding Window
# ----------------------------

class FastaParser:
    def __init__(self, sequences, motif_patterns):
        self.sequences = sequences
        self.motif_patterns = motif_patterns
        self.expanded_motifs = {}
        self.prepare_motifs()

    def prepare_motifs(self):
        self.expanded_motifs = {motif: expand_motif(motif) for motif in self.motif_patterns}

    def search(self):
        results = {}
        for seq_id, seq in self.sequences.items():
            results[seq_id] = {motif: [] for motif in self.motif_patterns}
            for motif in self.motif_patterns:
                expanded_set = self.expanded_motifs[motif]
                motif_length = len(motif)
                for i in range(len(seq) - motif_length + 1):
                    if seq[i:i + motif_length] in expanded_set:
                        results[seq_id][motif].append(i)
        return results
        

# ----------------------------
# Main Script Execution
# ----------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A motif sequence visualizing tool')
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input FASTA file name')
    parser.add_argument('-m', '--motif', type=str, required=True, help='Input Motif text file name')
    
    args = parser.parse_args()

    sequences = parse_fasta(args.fasta)
    motifs = parse_motifs(args.motif)

    fasta_parser = FastaParser(sequences, motifs)
    search_results = fasta_parser.search()
    print(search_results)
