#!/usr/bin/env python

import argparse


# Parse the args

parser = argparse.ArgumentParser(description = 'A motif sequence visualizing tool')
parser.add_argument('-f', '--fasta', type=str, required=True, help='Input FASTA file name')
parser.add_argument('-m', '--motif', type=str, required=True, help='Input Motif text file name')

args = parser.parse_args()

# Class initialization

class FastaParser:
    def __init__(self, ):
        '''
        Make a sliding window
        '''

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
    
    def expand_all():
        

        




        self.sequence
        # Make a dictionary where 
        motif_dict: dict = {}
        motif_set: set = {}
        motif_str: str = ""
        with open(args.motif, 'r') as m_fh:
            for i, line in enumerate(args.motif):
                line = line.strip()
                for char in line:
                    if char == 'Y':
                        # Make a copy of the previous motif_str
                        
                        # make a list or set or dict or something that holds all possible motifs for that line so I can makea  sliding window of that motif size to iterate over the file in the motif parser class.
                        for j in range(1):
                            motif_str += "C"
                            
                        
                        motif_set
                        if motif_dict[]
                        motif_dict


                    

