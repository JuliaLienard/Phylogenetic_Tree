# !/ usr/bin/ python3

'''
Usage: python SeqWithID.py fasta list.txt

'''

import sys
import re

fastaFile = sys.argv[1]
listContig = sys.argv[2]

ListToExclude = []

with open(listContig, "r") as IDlist:
	for line in IDlist:
		ID = line.strip()
		ListToExclude.append(ID)
print(f'Number of contig IDs is:',len(ListToExclude))


with open(fastaFile, "r") as fasta:
	header = ''
	HeaderWithContigID = 0 
	for line in fasta:
		line = line.strip()
		if line.startswith(">"):
			header = line
			if any(re.search(pattern, header) for pattern in ListToExclude):
				HeaderWithContigID += 1
print(f'Number of sequences with contig IDs to remove:',HeaderWithContigID)