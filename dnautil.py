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
	"this function will take a dictionary of sequences and report the length of them
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
	print("Longest sequence is %i basepairs %i in strand(s):%s"% (maxlen, len(maxnames), maxnames))
	print("Shortest sequence is %i basepairs %i in strand(s):%s"% (minlen, len(minnames), minnames))
	print("Checked across %i sequences"% (seqnum-1))
	return sequences

def lookforframe(dna):
	"Searches a string of basepairs to find a frame reference"
	dna = dna.upper()
	startcodon = "ATG"
	frameref = {1:False, 2:False, 3:False}
	for frame in range(0,3):
		for i in range(frame,len(dna),3):
			codon=dna[i:i+3]
			if codon == startcodon :
				frameref[frame+1] = True
				break
	validframes = 0
	totalvalids = 0
	if frameref[1] :
		validframe =1
		totalvalids +=1
	if frameref[2] :
		validframe =2
		totalvalids +=1
	if frameref[3] :
		validframe =3
		totalvalids +=1
	if totalvalids ==1:
		return validframe
	elif totalvalids==0:
		return "no valid frames"
	else:
		return "More than one valid frame"


def stripORFs(dna, frame=1):
	"searches a dna sequence for Open Reading Frames and returns a dictionary of the start position and the ORF"
	dna = dna.upper()
	startcodon = "ATG"
	stopcodons = ["TAA", "TAG", "TGA"]
	orfs = {}
	codonopened = False
	for i in range(frame-1, len(dna),3):
		codon=dna[i:i+3]
		if codon == startcodon :
			#if there was previously a startcodon with no stopcodon i'm ignoring it
			codonopened = True
			startbp = i+1
		elif codonopened and codon in stopcodons:
			codonopened = False #start looking for the next start
			orfs[startbp] = dna[startbp-1:i+3]
		else:
			pass
	return orfs

def fileORFs(fastafile):
	"takes a fastafile strips the sequences and then searches for ORFs"
	seqs = readfastafile(fastafile)
	fullseq={}
	for name,seq in seqs.items():
		frameref = lookforframe(seq)
		if type(frameref)="int":
			pass
		else:
			fullseq[name] = "No frame found"
