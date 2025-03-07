# Motif-Mark

Motif Mark is an object‑oriented Python tool designed to visualize motifs on biological sequences using pycairo, producing to-scale PNG images that clearly annotate sequence features such as exons, introns, and overlapping motifs. This project:

* Parses FASTA and motif files
* Dynamically resolves ambiguous bases
* Renders customizable visualizations with legends and color-coding

making it an ideal tool for researchers and bioinformatics enthusiasts. For detailed planning and design insights, please refer to my [OoCA Motif-Mark Plan](OoCA_Motif-mark.md).

The code: [Motif-Mark](motif-mark-oop.py)

The image: [Figure 1](Figure_1.svg)


## Usage
Command line options:

**-f (--fasta):** enter the FASTA file which contains all sequences you want to find motifs in. Each discrete sequence should be initialized with '>' inline header.

**-m (--motif):** enter the motif text file which contains all motif sequences you want to search for.
  One motif sequence per line.

  
