# !/ usr/bin/ python3

'''
Date: 2024-02-27
Author: Julia Lienard

Description: the program uses a list of contig IDs that need to be removed from a fasta file, 
including the DNA sequence.

Usage: python ExcudeContig.py fasta.fna (or.faa) ListContig.txt OutputFasta

'''

import sys
import re


fastaFile = sys.argv[1]
listContig = sys.argv[2]
fastaClean = sys.argv[3]


ListToExclude = []

with open(listContig, "r") as IDlist:
	for line in IDlist:
		ID = line.strip()
		ListToExclude.append(ID)



with open(fastaFile, "r") as fasta, open(fastaClean, "w") as output:
	header = ''
	sequence = ''
	for line in fasta:
		line = line.strip()
		if line.startswith(">"):
			header = line
			sequence = next(fasta)
			if not any(re.search(pattern, header) for pattern in ListToExclude):
				output.write(f'{header}\n{sequence}\n')
				sequence = ''
				header = ''
			


