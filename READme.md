# READme

The project was to do a phylogenetic analysis of Haemoproteus tartakovskyi (Ht) with other close relative parasite genomes (Plasmodium species).
The starting point was an already assembled genome of Haemoproteus tartakovskyi, contaminated with DNA of bird origin and: 
	1) Filter the contigs by size and GC content
	2a) Perform Gene prediction (using Genemark) for close related species (Plasmodium sp.)
	2b) Perform Gene prediction of the Ht genome (using Genemark) and then run BLAST to identify those of bird origin
	3) Clean the Ht genome to remove contigs containing predicted bird genes
	4) Run Gene prediction (using Genemark) on the cleaned Ht genome
	5) Compare the genome of Ht with close relatives (Plasmodium species)
	and the output group Toxoplasma gondii for size, gene numbers, GC content
	6) Identify orthologs using BUSCO between the Haemoproteus tartakovskyi genome and close relatives (Plasmodium species)
	and the output group Toxoplasma gondii
	7) Retrieve protein sequences of the orthologs to perform phylogenetic analyses
	8) Perform an alignment of the protein sequences for each gene orthologs, using clustalo
	9) Built a phylogenetic tree using raxmlHPC on each orthologs separately
	10) Make a consensus tree from all built trees using the phylip package


## STEP - 1: Filter the contigs by size and GC content

### get the data
```sh
# in Malaria directory
mkdir 1_rawdata
# copy plasmodiumGenomes and Haemoproteus tartakovskyi genome from course server
cd 1_rawdata
cp /home2/resources/binp29/Data/malaria/plasmodiumGenomes.tgz ./
cp /home2/resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz ./
# uncompress
tar -xvzf plasmodiumGenomes.tgz
gunzip Haemoproteus_tartakovskyi.raw.genome.gz
```

### Use removeScaffold.py to remove scaffolds < 3000 nucleotides and with GC% < in the Ht
```sh
cp /home2/resources/binp29/Data/malaria/removeScaffold.py ./
chmod +x removeScaffold.py  #make the script executable
# usage is removeScaffolds.py Haemoproteus_tartakovskyi.genome 35 Ht.genome 3000
# Based on the graph (in the presentation of the project) showing the GC content of the reads,\
# I choose 30% for the highest GC content that should be included. I also tried 35%
cd 1_rawdata
mkdir CleanHtgenome1
mv removeScaffolds.py ./CleanHtgenome1
cd CleanHtgenome1
./removeScaffold.py ../Haemoproteus_tartakovskyi.raw.genome 30 ./HtClean.genome 3000
./removeScaffold.py ../Haemoproteus_tartakovskyi.raw.genome 35 ./HtClean35.genome 3000
# number of scaffolds in the raw genome of Haemoproteus_tartakovskyi:
cat Haemoproteus_tartakovskyi.raw.genome | grep ^\> | wc -l # I get 15048
# number of scaffolds in the clean genomes of Haemoproteus tartakovskyi:
cat HtClean.genome | grep ^\> | wc -l # I get 2222
cat HtClean35.genome | grep ^\> | wc -l # 2343
# Percentage scaffolds left:
$ echo "2222*100/15048" | bc # 14% for HtClean.genome => I continue with this one
echo "2343*100/15048" | bc # 15% for HtClean35.genome
```

## STEP - 2a: Gene prediction for close related species (Plasmodium sp.)

```sh
# in Malaria directory
mkdir 2_GenePrediction
cd 2_GenePrediction
```
### Perform gene prediction for 1 genome using Genemark
from course tutorial : "To run GeneMark (gm) a license key is needed."
```sh
#key was downloaded from previous course, check:
ls ~/.gm_key
# Genemark is run using: gmes_petap.pl
nohup gmes_petap.pl --ES --sequence ../1_rawdata/Plasmodium_vivax.genome
# Genemark is run here using the self-training algorithm (--ES)
# The prediction files from all Plasmodium genomes are uploaded here:
/tmp
```

### Get the gene prediction file for Toxoplasma gondii from the course directory
```sh
# in 2_GenePrediction directory:
cp /home2/resources/binp29/Data/malaria/Tg.gff.gz ./
gunzip Tg.gff.gz # uncompress file
```

## STEP - 2b: Gene prediction of the Ht genome (using Genemark) and then run BLAST to identify those of bird origin
### Gene prediction for the clean Haemoproteus tartakovskyi genome (HtClean.genome)
```sh
cd 2_GenePrediction
mdkir 2_Ht
cd 2_Ht
nohup gmes_petap.pl --ES --min_contig 3000 --sequence ../../1_rawdata/CleanHtgenome/HtClean.genome
# the minimum size of the contigs is by default 50000 for humans, so here we need to specify --min_contig 3000 \
# otherwise we will get this error: "error, input sequence size is too small data/training.fna: 64494"
# rename the gtf file obtained:
mv genemark.gtf genemark.Ht.gtf
```

### get Fasta sequences from the gff file of HtClean.genome
Use gffParse.pl script from course server:
```sh
# in 2_GenePrediction/2_Ht:
cp /home2/resources/binp28/Data/gffParse.pl ./
chmod +x gffParse.pl # make the script executable
./gffParse.pl -h # get help
# usage : gffParse.pl -i fasta_file -g gff_file
./gffParse.pl -i ../../1_rawdata/CleanHtgenome1/HtClean.genome -g genemark.Ht.gtf
```

# header HtClean.genome:
>contig00001 GC=0.28 Length=64494
# first column genemark.Ht.gtf:
contig00001 GC=0.28 Length=64494
# Fix problem of compatible files (FASTA/gtf) to be able to run gffParse.pl:
cat genemark.Ht.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Ht2.gff

./gffParse.pl -i ../../1_rawdata/CleanHtgenome1/HtClean.genome -g genemark.Ht2.gff
# output files are gffParse.fna (fasta file containing the genes) and gffParse.log (log file)
## you can use a -p flag to also generate a .faa file!

## STEP - 3: Clean the Ht genome to remove contigs containing predicted bird genes
### Identify scaffolds with avian origin, using BLAST
```sh
echo $BLASTDB # to check if we have already the BLAST databases
```
```sh
# in /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/2_Ht/:
mkdir BLAST
cd BLAST
update_blastdb.pl swissprot # files downloaded: swissprot.tar.gz  swissprot.tar.gz.md5
# extract the swissprot database:
tar -xvzf swissprot.tar.gz

# Blastp using swissprot: I need to translate gffParse.fna into gffParse.faa (using the gffParse.pl, otherwise I used blastx
blastx -query gffParse.fna -db swissprot -evalue 1e-10 -out Ht.blastx -num_threads 10
```
### Identify the species related to birds:
```sh
# creating soft links for the taxomomy files from NBCI and Swissprot, respectively, into my the\
# following directory: /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/2_Ht/BLAST
ln -s /home2/resources/binp29/Data/malaria/taxonomy.dat
ln -s /home2/resources/binp29/Data/malaria/uniprot_sprot.dat

# Getting the datParser.py from course directory to identify Blast hits that are of bird origin
## setting a variable for the path to the fna
Htfna=/home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome1/HtClean.genome
python datParser.py Ht.blastx ../gffParse.fna taxonomy.dat uniprot_sprot.dat > Birdscaffolds.txt
cat Birdscaffolds.txt | wc -l
22 # scaffolds of bird origin that need to be removed from gffParse.fna
```
### Remove contigs with genes of bird origin in the Ht genome:
Python script ExcludeContig.py, to use the contig ID from Birdscaffolds.txt, parse the gffParse.fna, and print out only the headers and
corresponding sequence in the ouput file that do not belong to the contigs specified in Birdscaffolds.
```sh
# in /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/:
mkdir CleanHtgenome2
cd CleanHtgenome2
python /home/inf-27-2023/BINP29/3_Malaria/scripts/ExcludeContig.py \
/home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome1/HtClean.genome \
/home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/2_Ht/BLAST/Birdscaffolds.txt \
/home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome2/HtClean2.genome

# check it worked
cat /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/2_Ht/BLAST/Birdscaffolds.txt | wc -l # 22 contig ID to remove
cat /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome1/HtClean.genome | grep ^\> | wc -l # 2222 sequences before cleaning
cat /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome2/HtClean2.genome | grep ^\> | wc -l # 2200 sequences after cleaning
```

## STEP - 4: Gene prediction for the new clean Haemoproteus tartakovskyi genome (HtClean2.genome)
```sh
cd 2_GenePrediction
mdkir 2_HtClean
cd 2_HtClean
nohup gmes_petap.pl --ES --min_contig 3000 --sequence /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome2/HtClean2.genome &
# we specify again --min_contig 3000 \
# rename the gtf file obtained:
mv genemark.gtf genemark.HtClean2.gtf
# fix gff file for format (to be able to use gffParse.pl)
cat genemark.HtClean2.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Ht2.gff
```

## STEP - 5: Compare the genomes
### Create a table with information about the genomes (Genome size, Genes, Genomic GC)
```sh
# genome size:
cat HtClean2.genome | grep -v ^\> | tr -d "\n" | wc -c # 16563320
# GC content
cat HtClean2.genome | grep -v ^\> | tr -cd "cgCG" | wc -c # 4251511
echo "4251511*100/16563320" | bc # 25
# Genes number
cat genemark.HtClean2.gtf | cut -f 9 | awk -F ";" '{print$1}' | sort | uniq | wc -l # 3633
```
Do the same for all:
Species;genome size[bp];%GC;genes
Plasmodium_berghei;17954629;#23;7235
Plasmodium_cynomolgi;26181343;39;5787
Plasmodium_faciparum;23270305;19;5207
Plasmodium_knowlesi;23462346;37;4953
Plasmodium_vivax;27007701;42;5682
Plasmodium_yoelii;22222369;20;4889
Toxoplasma_gondii;128105889;80;15892

and collect in table

```sh
pip3 install pandas
# create genomesTable.py to make the table
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
```
```sh
chmod +x genomesTable.py
python genomesTable.py
```
#                    Species      Host Genome size (Mb) GC %  Genes
0  1         Plasmodium berghei   rodents               18   23   7235
1  2       Plasmodium cynomolgi  macaques               26   39   5787
2  3      Plasmodium falciparum    humans               23   19   5787
3  4        Plasmodium knowlesi    lemurs               23   37   5207
4  5           Plasmodium vivax    humans               27   42   4953
5  6          Plasmodium yoelii   rodents               22   20   5682
6  7  Haemoproteus tartakovskyi     birds               16   25   4889
7  8          Toxoplasma gondii    humans               80   52  15892


## STEP 6 - Identify orthologs
This can be done with different tools. Here use of proteinortho and BUSCO

### Print the protein sequences in fasta format for each genome, using gffParse.pl
#### Step1: Correct all gff files to use them with gffParser.pl
```sh
cat genemark.Pb.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Pb2.gff
cat genemark.Pc.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Pc2.gff
cat genemark.Pf.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Pf2.gff
cat genemark.Pk.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Pk2.gff
cat genemark.Pv.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Pv2.gff
cat genemark.Py.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > genemark.Py2.gff
```

#### Step2:
In Project folder
```sh
mkdir 3_IdentifOrtholgs
cd 3_IdentifOrtholgs
# set soft links for script and data needed for Ht genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/scripts/gffParse.pl
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/CleanHtgenome2/HtClean2.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/2_HtClean/genemark.Ht2.gff
perl gffParse.pl -b Ht -c -p -i HtClean2.genome -g genemark.Ht2.gff 
# set soft links for genome data needed for the other genomes
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_berghei.genome 
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_cynomolgi.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_faciparum.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_knowlesi.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_vivax.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Plasmodium_yoelii.genome
ln -s /home/inf-27-2023/BINP29/3_Malaria/1_rawdata/Toxoplasma_gondii.genome
# set soft links for gff files needed for the other genomes
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Pb2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Pc2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Pf2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Pk2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Pv2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/3_plasmodiumGenemark/genemark.Py2.gff
ln -s /home/inf-27-2023/BINP29/3_Malaria/2_GenePrediction/Tg.gff
# creating FASTA files for all other genomes:
perl gffParse.pl -b Pb -c -p -i Plasmodium_berghei.genome -g genemark.Pb2.gff
perl gffParse.pl -b Pc -c -p -i Plasmodium_cynomolgi.genome -g genemark.Pc2.gff
perl gffParse.pl -b Pf -c -p -i Plasmodium_faciparum.genome -g genemark.Pf2.gff
perl gffParse.pl -b Pk -c -p -i Plasmodium_knowlesi.genome -g genemark.Pk2.gff
perl gffParse.pl -b Pv -c -p -i Plasmodium_vivax.genome -g genemark.Pv2.gff
perl gffParse.pl -b Py -c -p -i Plasmodium_yoelii.genome -g genemark.Py2.gff
perl gffParse.pl -b Tg -c -p -i Toxoplasma_gondii.genome -g Tg.gff
# in the 3_IdentifOrtholgs directory:
mkdir FASTA
mv *.f*a ./FASTA/ # move all .faa and .fna in the FASTA folder
```

### Install software proteinortho (https://gitlab.com/paulklemm_PHD/proteinortho#installation)
in 3_IdentifOrtholgs directory:
```sh
conda create -n ProtOrtho # create conda environment
conda activate ProtOrtho
# then install proteinortho tool:
conda install -c bioconda proteinortho
# in the 3_IdentifOrtholgs/FASTA/ directory:
nohup proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &
```

### Use of BUSCO https://busco.ezlab.org/
Installation using conda (already done in previsous exercises)
```sh
conda activate QualBusco (BUSCO was installed in this environment)
# flag: -l LINEAGE, --lineage_dataset LINEAGE
# Specify the name of the BUSCO lineage to be used.
# - apicomplexa lineage is suggested
# in /home/inf-27-2023/BINP29/3_Malaria/3_IdentifOrthologs/FASTA
mkdir BUSCO
cd BUSCO
busco -i ../Pb.faa -o Pb -m prot -l apicomplexa
# I run for all of them: Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg uing a loop:
for sequences in *.faa;
do
	file=$(basename $sequences .faa);
	busco -i ${file}.faa -o ${file} -m prot -l apicomplexa
done
```

QUESTION-7 Compare how many BUSCOs (orthologues proteins) found in each proteome. Do the investigated parasites have close
to complete numbers of BUSCOs? 
=> I created a table BuscoSummary.txt with all numbers of BUSCO hits
Tg has only 3 complete single copy BUSCO hits, and 377 Complete/Duplicated, that we will have to use.

### Find the common BUSCO hits between the 8 species:
First Extract the BUSCO Ids for all Complete Single copy hits for all 7 species:
full_table.tsv found in the output directory of the BUSCO analysis for each species, called:
run_apicomplexa_odb10/
```sh
cat ./Ht/run_apicomplexa_odb10/full_table.tsv | grep "Complete" | awk -F "\t" '{print$1}' > BuscoHits_Ht.txt
# check the expected number of Complete Single copy hits (278) is found:
cat BuscoHits_Ht.txt | wc -l
# Extract the BUSCO Ids for all Complete Single copy hits in other species, the same way
# For Tg, we use the duplicated also:
cat ./Tg/run_apicomplexa_odb10/full_table.tsv | grep "Duplicated\|Complete" | awk -F "\t" '{print$1}'| sort | uniq > BuscoHits_Tg.txt
```

### Retrieve common BUSCO ID between the 8 species
Use the developped script to output the common elements of several lists (CommonElements.py)
in the script directory), by using a soft link in the \
```sh
/home/inf-27-2023/BINP29/3_Malaria/3_IdentifOrthologs/FASTA/BUSCO
ln -s ../../../../scripts/CommonElements.py 
## run the script on all BuscoHits_*.txt lists:
python CommonElements.py BuscoHits_*.txt BuscoCommonHit.txt
## the ouput is BuscoCommonHit.txt, which has 153 common BUSCO IDs between the 8 species
```
If Toxo is excluded, a higher number of common BUSCO (174) is found as this species is
more distantly related than the other parasites genomes.


## STEP 7: Retrieve protein fasta files from the 8 species for all 153 common BUSCO hits, in separate files (total 153):
### Step1: I rename the full_table.tsv obtained by the Ht BUSCO analysis, into Ht_full_table.tsv
Same is done for all species.
Then, I first create a table, using the script Busco2geneID in the script directory, that contains in 
column 1 the BuscoID (shared with all 8 species) and in the subsequent columns, the corresponding geneId for
each species.
```sh
In the home/inf-27-2023/BINP29/3_Malaria/3_IdentifOrthologs/FASTA/BUSCO:
mkdir BuscoTableAll
# I copy in this new directory all the "species"full_table.tsv obtained by the BUSCO analysis
cd BuscoTableAll
ln -s ln -s ../../../../scripts/Busco2geneID.py

python Busco2geneID.py BuscoCommonHit.txt Ht_full_table.tsv Pb_full_table.tsv Pc_full_table.tsv Pf_full_table.tsv Pk_full_table.tsv Pv_full_table.tsv Py_full_table.tsv Tg_full_table.tsv BuscoGeneID_All.txt
```

### Step2:
Now using the ouput BuscoGeneID_All.txt, the protein fasta sequences corresponding to each geneId in this table
need to be retrieve for the 8 species in a single file for each BUSCO Id.
Use of GeneIDtoSeq8Fasta.py in the script directory
```sh
# set soft link to this python script in /home/inf-27-2023/BINP29/3_Malaria/3_IdentifOrthologs/FASTA/:
ln -s /home/inf-27-2023/BINP29/3_Malaria/scripts/GeneIDtoSeq8Fasta.py
# set an output directory for the 153 generated files with the sequences:
mkdir output
python GeneIDtoSeq8Fasta.py BUSCO/BuscoTableAll/BuscoGeneID_All.txt Ht.faa Pb.faa Pc.faa Pf.faa Pk.faa Pv.faa Py.faa Tg.faa
ls -1 | wc -l # check we got 153 files
```

## STEP 8 - Alignment of the protein sequences using Clustalo

```sh
# Now we will align of 8 sequences in each of the 153 generated files.
# Install clustalo
conda create -n Alignment
conda activate Alignment
conda install -c bioconda clustalo raxml
```

```sh
# Run clustalo (I run it in the FASTA/output directory and then moved the aligned.faa files into the 4_alignments directory):
for sequences in *.faa;
do
	file=$(basename $sequences .faa);
	clustalo -i ${file}.faa -o ${file}_aligned.faa -v
done
```
### STEP 9 - Phylogenetic trees using raxmlHPC for each BUSCO
In the 4_alignments directory:
```sh
# create trees using raxmlHPC:

for alignements in *aligned.faa;
do
        file=$(basename "$alignements" _aligned.faa);
        raxmlHPC -s "${file}_aligned.faa" -n "${file}.tre" -o \
Tg -m PROTGAMMABLOSUM62 -p 12345
done

# move all output files *.tre and *.reduced in the output directory
```
### STEP 10 - Make consensus tree with all BUSCO-generated trees, using Phylip
```sh
#Install the phylip package to merge all individual trees:
conda install -c bioconda phylip

#I got different types of trees : bestTree, parsimonyTree and result files.
#I will merge all tree-files for each type into one.

#Create a directory for each type of tree, to make the consensus, using Phylip
mkdir ResultTree
mkdir bestTree
mkdir parsimonyTree

cat *result*.tre > ResultTree/Result.intree
cat *bestTree*.tre > bestTree/bestTree.intree
cat *parsimonyTree*.tre > parsimonyTree/parsimonyTree.intree

#Run Phylip in each corresponding directory:
phylip consense Result.intree
phylip consense bestTree.intree
phylip consense parsimonyTree.intree
```

```sh
#stdout:
Settings for this run:
 C         Consensus type (MRe, strict, MR, Ml):  Majority rule (extended)
 O                                Outgroup root:  No, use as outgroup species  1
 R                Trees to be treated as Rooted:  No
 T           Terminal type (IBM PC, ANSI, none):  ANSI
 1                Print out the sets of species:  Yes
 2         Print indications of progress of run:  Yes
 3                               Print out tree:  Yes
 4               Write out trees onto tree file:  Yes

Are these settings correct? (type Y or the letter for one to change)
O # I changed the outgroup to species 8 (Tg)

Consensus tree written to file "outtree"

Output written to file "outfile"

Done.
```

```sh
#Use of an online tool is possible to visualize the trees: https://itol.embl.de/upload.cgi
#by loading the info line in outfile
cat outfile
```
                                          +-------Pc
                                  +-100.0-|
                          +-124.0-|       +-------Pv
                          |       |
                  +--48.0-|       +---------------Pk
                  |       |
                  |       |               +-------Pb
          +--90.0-|       +---------149.0-|
          |       |                       +-------Py
  +-------|       |
  |       |       +-------------------------------Pf
  |       |
  |       +---------------------------------------Ht
  |
  +-----------------------------------------------Tg



```sh
cat outtree
```
((((((Pc:153.0,Pv:153.0):100.0,Pk:153.0):124.0,(Pb:153.0,Py:153.0):149.0):48.0,
Pf:153.0):90.0,Ht:153.0):153.0,Tg:153.0);

Here the tree obtained looks similar to the "true" species tree.
All Plasmodium species are clustering together, with the species falciparum on a separate branch than the others. P. falciparum has a very different GC content that the others.
The species having similar GC contents, seem to cluster together here.
Regarding the host range, it is hard to say on such a small sample size of species, but the species infecting rodents cluster together and the species infecting humans or monkeys also cluster together.








