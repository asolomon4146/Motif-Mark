#!/usr/bin/env python3

import argparse
import cairo

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
# Visualization with pycairo
# ----------------------------

def visualize_motifs_cairo(seq_id, sequence, search_results, output_file):
    """
    Create a visualization of the sequence with motifs marked using pycairo.
    - seq_id: identifier for the sequence.
    - sequence: nucleotide sequence (str).
    - search_results: dict mapping motif to a list of start indices.
    - output_file: PNG filename to save.
    """
    # Constants for drawing
    scale = 5  # pixels per base
    margin = 20
    seq_length = len(sequence)
    width = seq_length * scale + 2 * margin
    height = 150  # adjust as needed
    
    # Create a Cairo surface and context
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context(surface)
    
    # Fill background with white
    ctx.rectangle(0, 0, width, height)
    ctx.set_source_rgb(1, 1, 1)
    ctx.fill()
    
    # Draw the sequence as a horizontal line
    y_line = height // 2
    ctx.set_line_width(2)
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(margin, y_line)
    ctx.line_to(width - margin, y_line)
    ctx.stroke()
    
    # Define a list of colors (in RGB, 0-1 scale) for different motifs
    colors_list = [
        (1, 0, 0),     # red
        (0, 0, 1),     # blue
        (0, 0.8, 0),   # green
        (1, 0.5, 0),   # orange
        (0.5, 0, 0.5)  # purple
    ]
    
    # Draw motifs: for each motif, draw a rectangle at each occurrence
    for idx, (motif, positions) in enumerate(search_results.items()):
        color = colors_list[idx % len(colors_list)]
        for pos in positions:
            x = margin + pos * scale
            rect_width = len(motif) * scale
            # Vertical position: stagger motif markers above the line
            rect_y = y_line - 30 - idx * 15
            rect_height = 10
            ctx.rectangle(x, rect_y, rect_width, rect_height)
            ctx.set_source_rgb(*color)
            ctx.fill()
            
            # Optionally, annotate with the motif text
            ctx.set_source_rgb(0, 0, 0)
            ctx.set_font_size(10)
            ctx.move_to(x, rect_y - 2)
            ctx.show_text(motif)
    
    # Add a title to the image
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(14)
    ctx.move_to(margin, margin)
    ctx.show_text(f"Motif Occurrences in {seq_id}")
    
    # Save the image to a PNG file
    surface.write_to_png(output_file)

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

    # Output file prefix: same as FASTA file prefix
    output_prefix = args.fasta.rsplit('.', 1)[0]
    for seq_id, result in search_results.items():
        output_file = f"{output_prefix}_{seq_id}.png"
        visualize_motifs_cairo(seq_id, sequences[seq_id], result, output_file)
        print(f"Visualization saved to {output_file}")
