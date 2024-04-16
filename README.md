# OperonOrchestrator
This is code for handling sequences for the ribosomal operon database. Sequence input is from Mycocosm and NCBI.

Genomes in NCBI can be found: 
```
rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt .
```

Overview as of Oct. 16, 2023

```
##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.
## This is from column 25: 
  16538 archaea
1779617 bacteria
  16572 fungi
   6249 invertebrate
   7232 metagenomes
   7593 other (???)
   3604 plant
   1978 protozoa
   2854 vertebrate_mammalian
   3746 vertebrate_other
  80657 viral
```

TODO: Gather stats from Vsearch aligning on Saga for each genome
TODO: Select the reference sequences.
- should be longer than 4000 bp.
- If there are several copies they should be most similar to the others
