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
def printsequencelength(sequences):
	"this function will take a dictionary of sequences and report the length of them in a sorted manner"
	lengths={}
	for name,seq in sequences.items():
		lengths[name]=len(seq)
		print(name,len(seq))
	return lengths
def maxminseqs(sequences):
	"this function finds the max(s) & min(s) of a set of sequences and reports the names"
	maxlen = 0
	maxnames = []
	minlen = 0
	minnames = []
	seqnum = 1
	for name,seq in sequences.items():
		if seqnum ==1:
			#first record
			maxlen = len(seq)
			minlen = len(seq)
			maxnames=[name]
			minnames=[name]
		else:
			if len(seq) < minlen:
				#new shortest = clear old values
				minlen = len(seq)
				minnames=[name]
			elif len(seq) == minlen:
				#match shortest - add to list
				minnames.append(name)
			elif len(seq) > maxlen:
				#new longest = clear old values
				maxlen = len(seq)
				maxnames=[name]
			elif len(seq) == maxlen:
				#match longest - add to list
				maxnames.append(name)
		seqnum +=1
	print("Longest sequence is %i basepairs in strand(s):"% maxlen, maxnames)
	print("Shortest sequence is %i basepairs in strand(s):"% minlen, minnames)
	return sequences
