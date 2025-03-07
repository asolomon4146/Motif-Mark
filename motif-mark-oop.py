#!/usr/bin/env python

import argparse

# Helper function for file parsing

def parse_fasta(fasta_file):
    '''
    Parse FASTA file and return dict of sequences.
    Key: sequence identifier (minus the '>')
    Value: base sequence as a single string
    '''

    sequences: dict = {}
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
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = ''.join(seq_lines)
        
    return sequences

def parse_motifs(motif_file):
    '''
    Parse the motif file (one motif per line) and return a list of motifs
    '''

    motifs = []
    with open(motif_file, 'r') as motif_fh:
        for line in motif_fh:
            motif = line.strip()
            if motif:
                motifs.append(motif)
    return motifs


# Motif expansion class initialization


class MotifBranch:
    ambiguous_map: dict = {
        'Y' : ['T', 'C'], # What we really care about but I made it general for all ambiguous base calls
        'R' : ['A', 'G'],
        'S' : ['G', 'C'],
        'W' : ['A', 'T'],
        'N' : ['A', 'T', 'G', 'C']
    }
    
    def __init__(self, sequence = ''):
        '''
        A class to parse through the motif file and create individual motif objects
        '''
        self.sequence = sequence
        self.children = []

    def add_child(self, char):
        possibilities = MotifBranch.ambiguous_map.get(char, [char])
        for letter in possibilities: # create a new Motif branch object, iterate through T and C as a possibility
            child = MotifBranch(self.sequence + letter)
            self.children.append(child)

        return self.children
    
    def expand_all(self, motif_pattern):
        if not motif_pattern:
            return {self.sequence}
        
        char = motif_pattern[0]
        branches = self.add_child(char)
        results = set()

        for branch in branches:
            results.update(branch.expand_all(motif_pattern[1:]))
        return results
    
def expand_motif(motif):
    '''
    Given an ambiguous motif, return a set of all possible concrete motifs.
    '''
    root = MotifBranch()
    return root.expand_all(motif)

class FastaParser:
    def __init__(self, sequences, motif_patterns):
        '''
        sequences: dict of {seq_id: sequences}
        motif_patterns: list of ambiguous motif strings
        '''

        self.sequences = sequences
        self.motif_patterns = motif_patterns
        self.expanded_motifs = {}
        self.prepare_motifs()

    def prepare_motifs(self):
        '''
        expand each ambiguous motif into all possible motifs and store the results
        '''
        for motif in self.motif_patterns:
            expanded = expand_motif(motif)
            self.expanded_motifs[motif] = expanded

    def search(self):
        '''
        For each sequence and for each motif pattern, use a sliding window to find matches.
        Returns a dictionary of results: {seq_id: {motif: [list of start indices]}}
        '''
        search_dict: dict = {}
        for seq_id, seq in self.sequences.items():
            search_dict[seq_id] = {}
            for motif in self.motif_patterns:
                motif_length = len(motif)
                expanded_set = self.expanded_motifs[motif]
                search_dict[seq_id][motif] = []
                for i in range(len(seq) - motif_length + 1):
                    window = seq[i : i  + motif_length]
                    if window in expanded_set:
                        search_dict[seq_id][motif].append(i)


        return search_dict

if __name__ == '__main__':
    #-----------
    #argsparse
    #-----------
    parser = argparse.ArgumentParser(description = 'A motif sequence visualizing tool')
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input FASTA file name')
    parser.add_argument('-m', '--motif', type=str, required=True, help='Input Motif text file name')
    #-----------
    args = parser.parse_args()

    # root = MotifBranch()
    # complete_branches = root.expand_all("YYYY")
    # expanded_motifs = [branch.sequence for branch in complete_branches]
    # print(expanded_motifs)
    # print(len(expanded_motifs))
    # print(2**4)

    #parsing input files
    #-------------------
    sequences = parse_fasta(args.fasta)
    motifs = parse_motifs(args.motif)

    #Fastaparser instance
    fasta_parser = FastaParser(sequences, motifs)
    search_results_all = fasta_parser.search()
    print(search_results_all)
                    

