#!/usr/bin/env python3
import argparse
import cairo

# helper fn to parse fasta file
def parse_fasta(fasta_file):
    """parses a fasta file and returns a dict of sequences
    keeps original case so exons (upper) and introns (lower) can be told apart
    """
    seqs = {}
    current_id = None
    current_lines = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = "".join(current_lines)
                current_id = line[1:].split()[0]
                current_lines = []
            else:
                current_lines.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(current_lines)
    return seqs

# helper fn to parse motif file
def parse_motifs(motif_file):
    """reads a motif file (one motif per line) and returns a list of motifs in uppercase"""
    motifs = []
    with open(motif_file, 'r') as fh:
        for line in fh:
            motif = line.strip().upper()
            if motif:
                motifs.append(motif)
    return motifs

# class for expanding ambiguous motifs
class MotifBranch:
    # map ambiguous bases to possible letters
    ambig_map = {
        'Y': ['T', 'C'],
        'R': ['A', 'G'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'N': ['A', 'T', 'G', 'C']
    }
    
    def __init__(self, seq=""):
        self.seq = seq

    def expand_all(self, pattern):
        """recursively expand an ambiguous motif pattern into all possible sequences
        returns a set of sequences
        """
        if not pattern:
            return {self.seq}
        first = pattern[0]
        # get options for ambiguous base or default to itself
        opts = MotifBranch.ambig_map.get(first, [first])
        expanded_results = set()
        for letter in opts:
            new_seq = self.seq + letter
            expanded_results.update(MotifBranch(new_seq).expand_all(pattern[1:]))
        return expanded_results

def expand_motif(motif):
    return MotifBranch().expand_all(motif)

# parser class that searches for motifs in sequences
class FastaParser:
    def __init__(self, seqs, motifs):
        """seqs: dict of {id: sequence} and motifs is list of motif strings (upper)
        """
        self.seqs = seqs
        self.motifs = motifs
        self.expanded = {}
        self._prepare_motifs()

    def _prepare_motifs(self):
        """expand each motif to all possible unambiguous strings"""
        self.expanded = {motif: expand_motif(motif) for motif in self.motifs}

    def search(self):
        """search for motifs in every sequence (case-insensitive)
        returns dict {seq_id: {motif: [start_positions]}}
        """
        gene_motif_occurrences = {}
        for sid, seq in self.seqs.items():
            gene_motif_occurrences[sid] = {motif: [] for motif in self.motifs}
            seq_up = seq.upper()
            for motif in self.motifs:
                poss = self.expanded[motif]
                m_len = len(motif)
                for i in range(len(seq_up) - m_len + 1):
                    if seq_up[i:i+m_len] in poss:
                        gene_motif_occurrences[sid][motif].append(i)
        return gene_motif_occurrences

# visualizer class to draw the genes and motifs (output is SVG so its scalable)
class MotifVisualizer:
    def __init__(self, seqs, search_results, motifs, colors=None, scale=5, gene_height=100, gene_gap=20, margin=20):
        """
        seqs: dict of {id: sequence}
        search_results: dict {seq_id: {motif: [positions]}}
        motifs: list of original motif strings
        colors: list of RGB tuples (0-1) for motifs (if none, default colors used)
        scale: pixels per base
        gene_height: vertical space for gene drawing
        gene_gap: extra gap between genes to prevent overlaps
        margin: margin in pixels
        """
        self.seqs = seqs
        self.search_results = search_results
        self.motifs = motifs
        self.scale = scale
        self.gene_height = gene_height
        self.gene_gap = gene_gap
        self.margin = margin
        self.colors = colors if colors is not None else [
            (1, 0, 0),
            (0, 0, 1),
            (0, 0.8, 0),
            (1, 0.5, 0),
            (0.5, 0, 0.5)
        ]
        self.legend_height = 50
        self.max_len = max(len(seq) for seq in seqs.values())
        self.summary = {}

    def get_intervals(self, seq):
        """returns list of (start, end, region) for a seq
        region is 'exon' if letter is uppercase, 'intron' if lowercase
        """
        intervals = []
        if not seq:
            return intervals
        current = 'exon' if seq[0].isupper() else 'intron'
        start = 0
        for i, ch in enumerate(seq):
            region = 'exon' if ch.isupper() else 'intron'
            if region != current:
                intervals.append((start, i, current))
                start = i
                current = region
        intervals.append((start, len(seq), current))
        return intervals

    def assign_lanes(self, occs, m_len):
        """assign lanes to overlapping motif occurrances so they dont draw on top of each other
        returns dict {start: lane}
        """
        lanes = []
        lane_assign = {}
        for start in sorted(occs):
            placed = False
            for i, end in enumerate(lanes):
                if start >= end:
                    lanes[i] = start + m_len
                    lane_assign[start] = i
                    placed = True
                    break
            if not placed:
                lane_assign[start] = len(lanes)
                lanes.append(start + m_len)
        return lane_assign

    def classify(self, seq, start, m_len):
        """classify a motif occurrance as exon, intron or both
        """
        seg = seq[start:start+m_len]
        if all(ch.isupper() for ch in seg):
            return 'exon'
        elif all(ch.islower() for ch in seg):
            return 'intron'
        else:
            return 'both'

    def count_overlaps(self, positions):
        """count overlapping occurrances (if diff <= 1) in a sorted list
        """
        count = 0
        pos_sorted = sorted(positions)
        for i in range(len(pos_sorted)-1):
            if pos_sorted[i+1] - pos_sorted[i] <= 1:
                count += 1
        return count

    def draw(self, out_file):
        """draw the genes with motifs to an SVG file
        """
        width = self.max_len * self.scale + 2 * self.margin
        num_genes = len(self.seqs)
        spacing = self.gene_height + self.gene_gap
        height = num_genes * spacing + self.legend_height + self.margin

        surface = cairo.SVGSurface(out_file, width, height)
        ctx = cairo.Context(surface)

        # fill background white
        ctx.rectangle(0, 0, width, height)
        ctx.set_source_rgb(1, 1, 1)
        ctx.fill()

        self.summary = {}
        gene_idx = 0
        y_offset = self.margin

        for sid, seq in self.seqs.items():
            self.summary[sid] = {}
            base_y = y_offset + gene_idx * spacing + self.gene_height // 2

            # draw gene line (intron line) thicker
            ctx.set_line_width(3)
            ctx.set_source_rgb(0, 0, 0)
            start_x = self.margin
            end_x = self.margin + len(seq) * self.scale
            ctx.move_to(start_x, base_y)
            ctx.line_to(end_x, base_y)
            ctx.stroke()

            # draw exons (darker grey boxes)
            for (s, e, region) in self.get_intervals(seq):
                if region == 'exon':
                    exon_x = self.margin + s * self.scale
                    exon_w = (e - s) * self.scale
                    exon_h = 10
                    exon_y = base_y - exon_h // 2
                    ctx.rectangle(exon_x, exon_y, exon_w, exon_h)
                    ctx.set_source_rgb(0.5, 0.5, 0.5)
                    ctx.fill()

            # draw motifs and collect summary stats
            gene_occurrences = self.search_results.get(sid, {})
            for i, motif in enumerate(self.motifs):
                occs = gene_occurrences.get(motif, [])
                tot = len(occs)
                exon_ct = 0
                intron_ct = 0
                both_ct = 0
                for occ in occs:
                    typ = self.classify(seq, occ, len(motif))
                    if typ == 'exon':
                        exon_ct += 1
                    elif typ == 'intron':
                        intron_ct += 1
                    else:
                        both_ct += 1
                overlap_ct = self.count_overlaps(occs)
                self.summary[sid][motif] = {
                    'total': tot,
                    'exon': exon_ct,
                    'intron': intron_ct,
                    'both': both_ct,
                    'overlap': overlap_ct
                }
                lane_assign = self.assign_lanes(occs, len(motif))
                color = self.colors[i % len(self.colors)]
                for occ in occs:
                    lane = lane_assign[occ]
                    occ_start = occ
                    occ_end = occ + len(motif)
                    segments = []
                    curr_seg = None
                    seg_start = occ_start
                    for pos in range(occ_start, occ_end):
                        seg_type = 'exon' if seq[pos].isupper() else 'intron'
                        if curr_seg is None:
                            curr_seg = seg_type
                        elif seg_type != curr_seg:
                            segments.append((seg_start, pos, curr_seg))
                            seg_start = pos
                            curr_seg = seg_type
                    segments.append((seg_start, occ_end, curr_seg))
                    
                    motif_h = 8
                    # lane 0's bottom touches gene line
                    motif_y = base_y - motif_h - lane * (motif_h + 2)
                    for (seg_s, seg_e, seg_type) in segments:
                        x = self.margin + seg_s * self.scale
                        seg_w = (seg_e - seg_s) * self.scale
                        y_pos = motif_y
                        if seg_type == 'intron':
                            y_pos -= 2
                            ctx.set_source_rgba(*color, 0.6)
                        else:
                            ctx.set_source_rgba(*color, 1)
                        ctx.rectangle(x, y_pos, seg_w, motif_h)
                        ctx.fill()
                        ctx.set_line_width(0.5)
                        ctx.set_source_rgb(0, 0, 0)
                        ctx.rectangle(x, y_pos, seg_w, motif_h)
                        ctx.stroke()

            # label gene with bold font (size 25) moved 10px upward
            ctx.set_source_rgb(0, 0, 0)
            ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            ctx.set_font_size(25)
            # subtract 10 from y position
            gene_label_y = y_offset + gene_idx * spacing + 30 - 10
            ctx.move_to(self.margin, gene_label_y)
            ctx.show_text(sid)

            gene_idx += 1

        # calculate legend position so it fills the bottom white space and is 10px from border
        # using boxes 30x15 and key similar to motif legend
        legend_box_w = 30
        legend_box_h = 15
        legend_y = height - self.margin - 10 - (legend_box_h + 20)  # 20 extra for key area
        ctx.set_font_size(18)
        x_legend = self.margin
        for i, motif in enumerate(self.motifs):
            color = self.colors[i % len(self.colors)]
            ctx.rectangle(x_legend, legend_y, legend_box_w, legend_box_h)
            ctx.set_source_rgb(*color)
            ctx.fill()
            ctx.set_line_width(0.5)
            ctx.set_source_rgb(0, 0, 0)
            ctx.rectangle(x_legend, legend_y, legend_box_w, legend_box_h)
            ctx.stroke()
            ctx.move_to(x_legend + legend_box_w + 5, legend_y + legend_box_h)
            ctx.show_text(motif)
            x_legend += 150  # spacing between legend entries

        # draw key for gene structure (exon/intron)
        key_x = self.margin
        key_y = legend_y + legend_box_h + 5
        # Exon key
        ctx.rectangle(key_x, key_y, legend_box_w, legend_box_h)
        ctx.set_source_rgb(0.5, 0.5, 0.5)
        ctx.fill()
        ctx.set_source_rgb(0, 0, 0)
        ctx.rectangle(key_x, key_y, legend_box_w, legend_box_h)
        ctx.stroke()
        ctx.move_to(key_x + legend_box_w + 5, key_y + legend_box_h)
        ctx.show_text("Exon")
        key_x += 150
        # Intron key (line)
        ctx.set_line_width(3)
        ctx.move_to(key_x, key_y + legend_box_h // 2)
        ctx.line_to(key_x + legend_box_w, key_y + legend_box_h // 2)
        ctx.stroke()
        ctx.move_to(key_x + legend_box_w + 5, key_y + legend_box_h)
        ctx.show_text("Intron")

        surface.finish()
        print(f"visualization saved to {out_file}")

    def write_summary(self, summary_file):
        """writes summary stats to a tsv file with columns: Gene, Motif, Total, Exon, Intron, Both, Overlap
        """
        with open(summary_file, 'w') as fh:
            fh.write("Gene\tMotif\tTotal\tExon\tIntron\tBoth\tOverlap\n")
            for gene in self.summary:
                for motif in self.summary[gene]:
                    stats = self.summary[gene][motif]
                    line = f"{gene}\t{motif}\t{stats['total']}\t{stats['exon']}\t{stats['intron']}\t{stats['both']}\t{stats['overlap']}\n"
                    fh.write(line)
        print(f"summary saved to {summary_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='motif visualizer tool')
    parser.add_argument('-f', '--fasta', required=True, help='input fasta file')
    parser.add_argument('-m', '--motif', required=True, help='input motif file one motif per line')
    args = parser.parse_args()

    seqs = parse_fasta(args.fasta)
    motifs = parse_motifs(args.motif)

    parser_obj = FastaParser(seqs, motifs)
    gene_motif_occurrences = parser_obj.search()

    output_prefix = args.fasta.rsplit('.', 1)[0]
    out_svg = f"{output_prefix}.svg"
    summary_file = f"{output_prefix}_summary.tsv"

    visualizer = MotifVisualizer(seqs, gene_motif_occurrences, motifs)
    visualizer.draw(out_svg)
    visualizer.write_summary(summary_file)
