def repeatsearch(fastafile):
	"look for repeats of a certain length"
	"""
	A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) somewhere in the sequence. 
	Although repeats can occur on both the forward and reverse strands of the DNA sequence, 
	we will only consider repeats on the forward strand here. Also we will allow repeats to overlap themselves. 
	For example, the sequence ACACA contains two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3. 
	Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. 
	Your program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
	"""