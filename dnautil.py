#!/usr/bin/python
"""
Module of sequence manipulating functions
"""

def readfastafile(Fastafile):
	"this function will import a Fasta file and put the data into a dictionary"
	try:
		f = open(Fastafile)	
	except IOError:
		print("File %s does not exist!"% Fastafile)
	seqs={}
	for line in f:
		line = line.rstrip()
		if line.startswith('>'):
			words = line.split()
			name = words[0][1:]
			seqs[name]=''
		else:
			seqs[name] = seqs[name]+line
	f.close()
	return seqs
