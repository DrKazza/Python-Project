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
	"this function will take a dictionary of sequences and report the length of them"
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
		print("Found frame1")
		totalvalids +=1
	if frameref[2] :
		validframe =2
		print("Found frame2")
		totalvalids +=1
	if frameref[3] :
		validframe =3
		print("Found frame3")
		totalvalids +=1
	if totalvalids ==1:
		return validframe
	elif totalvalids==0:
		return  99
#		return "no valid frames"
	else:
		return  98
#		return "More than one valid frame"


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
		for i in range(1,4):
			subname = name + "Frame" + str(i)
			fullseq[subname] = stripORFs(seq,i)
	return fullseq

def longestORFinfile(fastafile,frame):
	"takes a fastafile strips the sequences and then searches for ORFs on a specific frame"
	seqs = readfastafile(fastafile)
	fullseq={}
	longestORF = 0
	seqname = []
	Longestitem = {}
	for name,seq in seqs.items():
		fullseq[name] = stripORFs(seq,frame)
	for name,orfs in fullseq.items():
		for startbp,bases in orfs.items():
			thislength = len(bases)
			if thislength > longestORF:
				seqname = [name]
				tmp = startbp
				longestORF = thislength
			elif thislength== longestORF:
				seqname.append(name)
			else:
				pass
	print(startbp)
	Longestitem[longestORF] = seqname
	return Longestitem

def longestORFinname(seqname, fastafile):
	"finds a sequence in a file then will find the longest ORF the start position and the frame which it's on"
	seqs = readfastafile(fastafile)
	strippedorfs = {}
	longestORF = []
	lengthlongest = 0
	longeststart = []
	longestframe = []
	for i in range(1,4):
		thisframe = stripORFs (seqs[seqname],i)
		for startbp, bases in thisframe.items():
			if len(bases) > lengthlongest:
				lengthlongest = len(bases)
				longestORF = [bases]
				longeststart = [startbp]
				longestframe = [i]
			elif len(bases) == lengthlongest:
				longestORF.append(bases)
				longeststart.append[startbp]
				longestframe.append[i]
			else:
				pass
	print("Longest ORF is %ibp starting at base %s on frame %s"% (lengthlongest, longeststart, longestframe))
	print(longestORF)
	return

def findrepeats(fastafile, replen):
	"Reads a fastafile and looks for all repeats, reports the number of repeats of length replen and the most frequent repeat"
	seqs = readfastafile(fastafile)
	repeatlist = {}
	maxrepeat = []
	maxrepno = 0
	for name,seq in seqs.items():
		for i in range(0,len(seq)-replen+1):
			searchstring = seq[i:replen+i]
			count = 0
			for seqname,seqbps in seqs.items():
				occur = seqbps.count(searchstring)
				count = count + occur
			if count > 1:
				repeatlist[searchstring] = count
				if count>maxrepno:
					maxrepno = count
					maxrepeat = [searchstring]
				elif count == maxrepno and searchstring not in maxrepeat:
					maxrepeat.append(searchstring)
				else:
					pass
	print("The most repeated sequence of %ibps occurs %i times and is %s"% (replen, maxrepno, maxrepeat))
	return repeatlist
