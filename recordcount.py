def recordcount(fastafile):
	"How many records are in the file?"
	
	"""
	A record in a FASTA file is defined as a single-line header, followed by lines of sequence data. 
	The header line is distinguished from the sequence data by a greater-than (>) symbol in the first column. 
	The word following the > symbol is the identifier of the sequence, 
	and the rest of the line is an optional description of the entry. There should be no space between the > and the first letter of the identifier.
	"""

	i = 0
	a = 0
