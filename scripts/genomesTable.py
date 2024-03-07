#!/ usr/bin/ python3

import pandas as pd

data = {
	"#": [1, 2, 3, 4, 5, 6, 7, 8],
	"Species": ["Plasmodium berghei", "Plasmodium cynomolgi", "Plasmodium falciparum", "Plasmodium knowlesi", "Plasmodium vivax", "Plasmodium yoelii", "Haemoproteus tartakovskyi", "Toxoplasma gondii"],
	"Host": ["rodents", "macaques", "humans", "lemurs", "humans", "rodents", "birds", "humans"],
	"Genome size (Mb)": ["18", "26", "23", "23", "27", "22", "16", "80"],
	"GC %": ["23", "39", "19", "37", "42", "20", "25", "52"],
	"Genes": ["7235", "5787", "5787", "5207", "4953", "5682", "4889", "15892"]
}

df = pd.DataFrame(data)

print(df)
