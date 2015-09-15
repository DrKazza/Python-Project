def readfastafile(Fastafile)
	"""
	reads a FastA file and builds a dictionary with the sequences and Names
	"""
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
return seqs
f.close()