# !/ usr/bin/ python3

'''
Date: 2024-02-29
Author: Julia Lienard

Description: parse a table with gene Id and recover from 8 fasta files, the corresponding header and sequence.

Usage: python GeneIDtoSeq8Fasta.py BuscoGeneID_All.txt Ht.faa Pb.faa Pc.faa Pf.faa Pk.faa Pv.faa Py.faa Tg.faa


'''

import sys

GeneIDList = sys.argv[1]
FastaHt = sys.argv[2]
FastaPb = sys.argv[3]
FastaPc = sys.argv[4]
FastaPf = sys.argv[5]
FastaPk = sys.argv[6]
FastaPv = sys.argv[7]
FastaPy = sys.argv[8]
FastaTg = sys.argv[9]

with open(GeneIDList, "r") as genelist:
	for BuscoId in genelist:
		col = BuscoId.strip().split("\t")
		filename = col[0] + ".faa"
		GeneId_Ht = col[1]
		GeneId_Pb = col[2]
		GeneId_Pc = col[3]
		GeneId_Pf = col[4]
		GeneId_Pk = col[5]
		GeneId_Pv = col[6]
		GeneId_Py = col[7]
		GeneId_Tg = col[8]
		SequenceList = {}
		with open(FastaHt, "r") as SeqHt:
			for line in SeqHt:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Ht:
						sequence = next(SeqHt)
						SequenceList["Ht\t"+ geneId] = sequence
						break
		with open(FastaPb, "r") as SeqPb:
			for line in SeqPb:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Pb:
						sequence = next(SeqPb)
						SequenceList["Pb\t"+ geneId] = sequence
						break
		with open(FastaPc, "r") as SeqPc:
			for line in SeqPc:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Pc:
						sequence = next(SeqPc)
						SequenceList["Pc\t"+ geneId] = sequence
						break
		with open(FastaPf, "r") as SeqPf:
			for line in SeqPf:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Pf:
						sequence = next(SeqPf)
						SequenceList["Pf\t"+ geneId] = sequence
						break
		with open(FastaPk, "r") as SeqPk:
			for line in SeqPk:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Pk:
						sequence = next(SeqPk)
						SequenceList["Pk\t"+ geneId] = sequence
						break

		with open(FastaPv, "r") as SeqPv:
			for line in SeqPv:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Pv:
						sequence = next(SeqPv)
						SequenceList["Pv\t"+ geneId] = sequence
						break
		with open(FastaPy, "r") as SeqPy:
			for line in SeqPy:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Py:
						sequence = next(SeqPy)
						SequenceList["Py\t"+ geneId] = sequence
						break
		with open(FastaTg, "r") as SeqTg:
			for line in SeqTg:
				if line.startswith(">"):
					header = line[1:].rstrip().split("\t")
					geneId = header[0].strip()
					if geneId == GeneId_Tg:
						sequence = next(SeqTg)
						SequenceList["Tg\t"+ geneId] = sequence
						break


		with open("output/" + filename, "w") as outputfile:
			for key, value in SequenceList.items():
				outputfile.write(f'>{key}\n{value}\n')
