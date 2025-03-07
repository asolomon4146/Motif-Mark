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
    NOTE: The original case is preserved so that exons (uppercase) and introns (lowercase) can be distinguished.
    """
    sequences = {}
    with open(fasta_file, 'r') as fa_fh:
        seq_id = None
        seq_lines = []
        for line in fa_fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq_lines)
    return sequences

def parse_motifs(motif_file):
    """
    Parse the motif file (one motif per line) and return a list of motifs.
    The motifs are converted to uppercase (for matching and ambiguous expansion).
    """
    motifs = []
    with open(motif_file, 'r') as motif_fh:
        for line in motif_fh:
            motif = line.strip().upper()
            if motif:
                motifs.append(motif)
    return motifs

# ----------------------------
# Motif Expansion Class (for ambiguous nucleotides)
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
        Returns a set of expanded motif strings.
        """
        if not motif_pattern:
            return {self.sequence}
        char = motif_pattern[0]
        letters = MotifBranch.ambiguous_map.get(char, [char])
        results = set()
        for letter in letters:
            new_seq = self.sequence + letter
            results.update(MotifBranch(new_seq).expand_all(motif_pattern[1:]))
        return results

def expand_motif(motif):
    return MotifBranch().expand_all(motif)

# ----------------------------
# FASTA Parser and Motif Searcher (using sliding window)
# ----------------------------

class FastaParser:
    def __init__(self, sequences, motif_patterns):
        """
        sequences: dictionary of {seq_id: sequence (with case preserved)}
        motif_patterns: list of motifs (in uppercase, may contain ambiguous bases)
        """
        self.sequences = sequences
        self.motif_patterns = motif_patterns
        self.expanded_motifs = {}
        self.prepare_motifs()

    def prepare_motifs(self):
        """Expand each ambiguous motif into all possible unambiguous sequences."""
        self.expanded_motifs = {motif: expand_motif(motif) for motif in self.motif_patterns}

    def search(self):
        """
        Search for each motif in every sequence.
        Returns a nested dictionary: {seq_id: {motif: [start_positions]}}.
        Matching is case-insensitive.
        """
        results = {}
        for seq_id, seq in self.sequences.items():
            results[seq_id] = {motif: [] for motif in self.motif_patterns}
            seq_upper = seq.upper()
            for motif in self.motif_patterns:
                expanded_set = self.expanded_motifs[motif]
                motif_len = len(motif)
                for i in range(len(seq_upper) - motif_len + 1):
                    if seq_upper[i:i + motif_len] in expanded_set:
                        results[seq_id][motif].append(i)
        return results

# ----------------------------
# Visualization with pycairo (SVG output)
# ----------------------------

class MotifVisualizer:
    def __init__(self, sequences, search_results, motif_patterns, colors=None, scale=5, gene_height=130, gene_gap=20, margin=17):
        """
        sequences: dict {seq_id: sequence}
        search_results: dict {seq_id: {motif: [start_positions]}}
        motif_patterns: list of original motif strings (with ambiguous bases)
        colors: optional list of RGB tuples (values 0-1) for motifs
        scale: pixels per base (width scaling factor)
        gene_height: vertical space allocated for gene drawing (affects exon/intron display)
        gene_gap: additional vertical gap between gene rows to prevent overlapping staggered motifs
        margin: margin in pixels around the drawing
        """
        self.sequences = sequences
        self.search_results = search_results
        self.motif_patterns = motif_patterns
        self.scale = scale
        self.gene_height = gene_height
        self.gene_gap = gene_gap
        self.margin = margin
        self.colors = colors if colors is not None else [
            (1, 0, 0),     # red
            (0, 0, 1),     # blue
            (0, 0.8, 0),   # green
            (1, 0.5, 0),   # orange
            (0.5, 0, 0.5)  # purple
        ]
        self.legend_height = 50  # this will be scaled up later in drawing
        self.max_seq_length = max(len(seq) for seq in sequences.values())
        # To store summary stats for output to file
        self.summary_stats = {}

    def get_exon_intron_intervals(self, seq):
        '''
        Determines exon and intron intervals.
        Exons are defined as contiguous uppercase letters,
        while introns are contiguous lowercase letters.
        Returns a list of tuples: (start, end, region_type)
        '''
        intervals: list = []
        if not seq:
            return intervals
        current_type = 'exon' if seq[0].isupper() else 'intron'
        start = 0
        for i, char in enumerate(seq):
            char_type = 'exon' if char.isupper() else 'intron'
            if char_type != current_type:
                intervals.append((start, i, current_type))
                start = i
                current_type = char_type
        intervals.append((start, len(seq), current_type))
        return intervals

    def assign_motif_lanes(self, occurrences, motif_length):
        """
        For a given motif (list of start positions) assign lane numbers so that
        overlapping occurrences are drawn in separate vertical "lanes".
        Returns a dictionary mapping occurrence start position to lane number.
        """
        lanes = []  # each lane holds the end position of the last occurrence placed in it
        lane_assignment = {}
        for start in sorted(occurrences):
            placed = False
            for lane_index, lane_end in enumerate(lanes):
                if start >= lane_end:  # no overlap with previous occurrence in this lane
                    lanes[lane_index] = start + motif_length
                    lane_assignment[start] = lane_index
                    placed = True
                    break
            if not placed:
                lane_index = len(lanes)
                lanes.append(start + motif_length)
                lane_assignment[start] = lane_index
        return lane_assignment

    def classify_occurrence(self, seq, start, motif_length):
        '''
        Classify a motif occurrence as 'exon', 'intron', or 'both'.
        '''
        segment = seq[start:start+motif_length]
        if all(c.isupper() for c in segment):
            return 'exon'
        elif all(c.islower() for c in segment):
            return 'intron'
        else:
            return 'both'

    def count_overlaps(self, positions):
        """
        Count overlapping occurrences in a sorted list of positions.
        Two occurrences are considered overlapping if the difference is <= 1.
        """
        overlaps = 0
        sorted_positions = sorted(positions)
        for i in range(len(sorted_positions)-1):
            if sorted_positions[i+1] - sorted_positions[i] <= 1:
                overlaps += 1
        return overlaps

    def draw(self, output_file):
        """
        Draw all genes and their motif occurrences on a single image.
        Gene structure (exons and introns) is drawn to scale.
        Overlapping motif occurrences are staggered, and if a motif occurrence spans
        both exon and intron regions, it is subdivided and drawn with slight vertical offsets and opacity differences.
        The drawing is output as an SVG file.
        """
        width = self.max_seq_length * self.scale + 2 * self.margin
        num_genes = len(self.sequences)
        gene_spacing = self.gene_height + self.gene_gap
        height = num_genes * gene_spacing + self.legend_height + self.margin

        # Create an SVG surface
        surface = cairo.SVGSurface(output_file, width, height)
        ctx = cairo.Context(surface)

        # Fill background white
        ctx.rectangle(0, 0, width, height)
        ctx.set_source_rgb(1, 1, 1)
        ctx.fill()

        gene_index = 0
        gene_y_offset = self.margin

        # Initialize summary stats dictionary
        self.summary_stats = {}

        for seq_id, seq in self.sequences.items():
            self.summary_stats[seq_id] = {}
            # Center position for the gene line in this gene row
            base_y = gene_y_offset + gene_index * gene_spacing + self.gene_height // 2

            # Draw gene line (representing introns) with increased thickness
            ctx.set_line_width(3)
            ctx.set_source_rgb(0, 0, 0)
            start_x = self.margin
            end_x = self.margin + len(seq) * self.scale
            ctx.move_to(start_x, base_y)
            ctx.line_to(end_x, base_y)
            ctx.stroke()

            # Draw exons as darker grey boxes over the gene line
            intervals = self.get_exon_intron_intervals(seq)
            for (start, end, region_type) in intervals:
                if region_type == 'exon':
                    exon_x = self.margin + start * self.scale
                    exon_width = (end - start) * self.scale
                    exon_height = 10  # fixed height for exon boxes
                    exon_y = base_y - exon_height // 2
                    ctx.rectangle(exon_x, exon_y, exon_width, exon_height)
                    ctx.set_source_rgb(0.5, 0.5, 0.5)  # darker grey for exons
                    ctx.fill()

            # Draw motif occurrences and collect summary stats
            gene_results = self.search_results.get(seq_id, {})
            for i, motif in enumerate(self.motif_patterns):
                occ_list = gene_results.get(motif, [])
                # Initialize summary stats for this motif
                total = len(occ_list)
                exon_count = 0
                intron_count = 0
                both_count = 0
                for occ in occ_list:
                    classification = self.classify_occurrence(seq, occ, len(motif))
                    if classification == 'exon':
                        exon_count += 1
                    elif classification == 'intron':
                        intron_count += 1
                    else:
                        both_count += 1
                overlapping = self.count_overlaps(occ_list)
                self.summary_stats[seq_id][motif] = {
                    'total': total,
                    'exon': exon_count,
                    'intron': intron_count,
                    'both': both_count,
                    'overlap': overlapping
                }
                lane_assignment = self.assign_motif_lanes(occ_list, len(motif))
                color = self.colors[i % len(self.colors)]
                for start in occ_list:
                    lane = lane_assignment[start]
                    occ_start = start
                    occ_end = start + len(motif)
                    # Break the occurrence into segments based on exon/intron
                    segments = []
                    current_seg_type = None
                    seg_start = occ_start
                    for pos in range(occ_start, occ_end):
                        char_type = 'exon' if seq[pos].isupper() else 'intron'
                        if current_seg_type is None:
                            current_seg_type = char_type
                        elif char_type != current_seg_type:
                            segments.append((seg_start, pos, current_seg_type))
                            seg_start = pos
                            current_seg_type = char_type
                    segments.append((seg_start, occ_end, current_seg_type))

                    # Motif dimensions
                    motif_height = 8
                    # Adjusted motif_base_y: lane 0's bottom touches the gene line.
                    motif_base_y = base_y - motif_height - lane * (motif_height + 2)

                    # Draw each segment of the motif occurrence
                    for (seg_start, seg_end, seg_type) in segments:
                        x = self.margin + seg_start * self.scale
                        seg_width = (seg_end - seg_start) * self.scale
                        y_offset = motif_base_y
                        if seg_type == 'intron':
                            y_offset -= 2  # shift intron portions slightly upward
                            ctx.set_source_rgba(*color, 0.6)
                        else:
                            ctx.set_source_rgba(*color, 1)
                        ctx.rectangle(x, y_offset, seg_width, motif_height)
                        ctx.fill()
                        # Draw a thinner outline around the motif segment
                        ctx.set_line_width(0.5)
                        ctx.set_source_rgb(0, 0, 0)
                        ctx.rectangle(x, y_offset, seg_width, motif_height)
                        ctx.stroke()

            # Label gene (e.g., gene ID) with a larger font
            ctx.set_source_rgb(0, 0, 0)
            ctx.set_font_size(16)  # Increased gene label font size
            label_y = gene_y_offset + gene_index * gene_spacing + 20
            ctx.move_to(self.margin, label_y)
            ctx.show_text(seq_id)

            gene_index += 1

        # Draw legend at the bottom (scaled up by 150%)
        legend_y = gene_y_offset + num_genes * gene_spacing + 10
        ctx.set_font_size(18)  # increased legend font size (150% of 12)
        x_legend = self.margin

        # Motif legend: show the original ambiguous motif with its color
        for i, motif in enumerate(self.motif_patterns):
            color = self.colors[i % len(self.colors)]
            ctx.rectangle(x_legend, legend_y, 30, 15)  # larger colored box
            ctx.set_source_rgb(*color)
            ctx.fill()
            ctx.set_line_width(0.5)
            ctx.set_source_rgb(0, 0, 0)
            ctx.rectangle(x_legend, legend_y, 30, 15)
            ctx.stroke()
            ctx.move_to(x_legend + 35, legend_y + 15)
            ctx.show_text(motif)
            x_legend += 140  # spacing between legend entries

        # Key for gene structure
        key_x = self.margin
        key_y = legend_y + 25
        # Exon key (darker grey box)
        ctx.rectangle(key_x, key_y, 30, 15)
        ctx.set_source_rgb(0.5, 0.5, 0.5)
        ctx.fill()
        ctx.set_source_rgb(0, 0, 0)
        ctx.rectangle(key_x, key_y, 30, 15)
        ctx.stroke()
        ctx.move_to(key_x + 35, key_y + 15)
        ctx.show_text("Exon")
        key_x += 140
        # Intron key (thicker line)
        ctx.set_line_width(3)
        ctx.move_to(key_x, key_y + 7)
        ctx.line_to(key_x + 30, key_y + 7)
        ctx.stroke()
        ctx.move_to(key_x + 35, key_y + 15)
        ctx.show_text("Intron")

        # Finish the SVG surface
        surface.finish()
        print(f"Visualization saved to {output_file}")

    def write_summary_stats(self, summary_file):
        """
        Write the summary statistics (per gene per motif) to a TSV file.
        Columns: Gene, Motif, Total, Exon, Intron, Both, Overlap
        """
        with open(summary_file, 'w') as outfh:
            header = "\t".join(["Gene", "Motif", "Total", "Exon", "Intron", "Both", "Overlap"])
            outfh.write(header + "\n")
            for gene in self.summary_stats:
                for motif in self.summary_stats[gene]:
                    stats = self.summary_stats[gene][motif]
                    line = "\t".join([gene, motif,
                                      str(stats['total']),
                                      str(stats['exon']),
                                      str(stats['intron']),
                                      str(stats['both']),
                                      str(stats['overlap'])])
                    outfh.write(line + "\n")
        print(f"Summary statistics saved to {summary_file}")

# ----------------------------
# Main Functionality
# ----------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A motif sequence visualizing tool')
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-m', '--motif', type=str, required=True, help='Input motif file (one motif per line)')
    args = parser.parse_args()

    # Parse input files
    sequences = parse_fasta(args.fasta)
    motifs = parse_motifs(args.motif)

    # Create FastaParser instance and search for motifs
    fasta_parser = FastaParser(sequences, motifs)
    search_results = fasta_parser.search()

    # Output file: use the FASTA file prefix (e.g. Figure_1.fasta -> Figure_1.svg)
    output_prefix = args.fasta.rsplit('.', 1)[0]
    output_file = f"{output_prefix}.svg"
    summary_file = f"{output_prefix}_summary.tsv"

    # Create visualizer instance, draw the image, and write summary stats to a file
    visualizer = MotifVisualizer(sequences, search_results, motifs)
    visualizer.draw(output_file)
    visualizer.write_summary_stats(summary_file)
