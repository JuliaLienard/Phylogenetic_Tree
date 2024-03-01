# !/ usr/bin/ python3

'''
Date: 2024-02-29
Author: Julia Lienard
Create a tab-separated able with BUSCO Id; corrresponding geneId species 1: geneID species2...


Usage:
 python BuscoTogeneIDfinal.py ID_Table BuscoTables BuscoTablespecies2 BuscoTablespecies3 BuscoTablespecies4 \
 BuscoTablespecies5 BuscoTablespecies6 BuscoTablespecies7 BuscoTablespecies8 output.txt
'''

import sys

BuscoIDList = sys.argv[1]
TableHt = sys.argv[2]
TablePb = sys.argv[3]
TablePc = sys.argv[4]
TablePf = sys.argv[5]
TablePk = sys.argv[6]
TablePv = sys.argv[7]
TablePy = sys.argv[8]
TableTg = sys.argv[9]
OutputTable = sys.argv[10]

# all the input files should be in the input folder

with open(BuscoIDList, "r") as BuscoID, open(OutputTable, "w") as output:
	output.write(f'BuscoID\tHt_geneID\tPb_geneID\tPc_geneID\tPf_geneID\tPk_geneID\tPv_geneID\tPy_geneID\tTg_geneID\n')
	for lineID in BuscoID:
		ID = lineID.strip()
		GeneIdList = []
		with open(TableHt, "r") as Ht:
			for rowHT in Ht:
				if not rowHT.startswith("#"):
					col = rowHT.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break # to avoid scanning the whole list even after finding the right element
		with open(TablePb, "r") as Pb:
			for rowPb in Pb:
				if not rowPb.startswith("#"):
					col = rowPb.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break # to avoid scanning the whole list even after finding the right element

		with open(TablePc, "r") as Pc:
			for rowPc in Pc:
				if not rowPc.startswith("#"):
					col = rowPc.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break 
		with open(TablePf, "r") as Pf:
			for rowPf in Pf:
				if not rowPf.startswith("#"):
					col = rowPf.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break
		with open(TablePk, "r") as Pk:
			for rowPk in Pk:
				if not rowPk.startswith("#"):
					col = rowPk.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break
		with open(TablePv, "r") as Pv:
			for rowPv in Pv:
				if not rowPv.startswith("#"):
					col = rowPv.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break
		with open(TablePy, "r") as Py:
			for rowPy in Py:
				if not rowPy.startswith("#"):
					col = rowPy.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break
		with open(TableTg, "r") as Tg:
			for rowTg in Tg:
				if not rowTg.startswith("#"):
					col = rowTg.strip().split("\t")
					buscoID = col[0].strip()
					if len(col) > 2:
						GeneId = col[2].strip()
						if ID == buscoID:
							GeneIdList.append(GeneId)
							break

			StringList = "\t".join(GeneIdList)
			output.write(f'{ID}\t{StringList}\n')

