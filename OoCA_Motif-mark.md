
## Classes
1. FASTA sequence parser
	1. Creates 1 or more `sequence` objects which could have the following attributes:
		1. Sequence ID/Name:
		2. Nucleotide String: The full sequence.
		3. Annotation Information:
			1. Exon/Intron flags
			2. Start/End pos for each exon/intron
	2. With multiple sequence objects their positions could be read by the motif parser to find motifs from the motif file within the sequences
	3. Methods:
		1. Parse:
			1. Reads the FASTA file
			2. extracts headers and sequences
			3. Returns a list of sequence objects
		2. Annotate:
			1. Assign whether a region is an exon or intron
			2. Store the start/end pos of these features
	4. 
2. Motif parser
	1. Creates 1 or more `motif` objects which could have the following attributes:
		1. The motif string
		2. properties (handling ambiguous bases)
	2. Could have the following methods:
		1. matching the motif in a sequence object (see above in FASTA sequence parser)
			1. Could return the start/end position 
		2. Resolve ambiguity: convert ambiguous bases to a regex to match all potential nucleotide combinations
3. Visualizer
	1. creates visualizations as rectangle or line objects using cairo with the following possible attributes:
		1. X and y start euclidean start position and their lengths/widths
		2. stroke or fill
		3. name of rectangle
		4. color
		5. Scale factor
	2. Methods:
		1. Checks with the motif and FASTA parsers to decide what color to make the rectangles depending on which motif is used. Could make it possibly generate hex codes dynamically to accommodate unlimited potential motifs but might not be necessary (possibly necessary for the `YYYYYYYYYY` motif)
		2. Canvas initialize
			1. initialize the pycairo canvas with dimensions
		3. Draw sequence
			1. Draw the base line representing the full sequence
			2. Annotate exon/intron boundaries
		4. Draw motif
			1. Draw a rectangle at the corresponding positions
		5. Create a legend
			1. Include:
				1. Color
				2. Annotation start/end pos of exon/intron regions
		6. Save file
			1. Save as a png with correct naming convention
