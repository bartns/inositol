*Tested on Ubuntu 18.04*

**Prerequisites:**
- Tools: curl, parallel
- Python3
- Python3 Packages: numpy, pandas, matplotlib, seaborn
- DIAMOND v0.9.22: http://www.diamondsearch.org
- KRAKEN v2.0.9-beta: https://ccb.jhu.edu/software/kraken2/
- BRACKEN v2.5: https://ccb.jhu.edu/software/bracken/
- The files in this folder


 File explanations:

- **RAST_mgmIDs:** RAST data files used for this analyses
- **SRR_datasets:** corresponding SRR dataset IDs from RAST
- **rename_RastID_to_SRR.sh:** Rename file with RASTid to it's corresponding SRR id (only for this dataset)
- **totalReads_samples.tsv:** Total reads in sample. (RAST screened and passed fastq)
  

- **Anaerostipes_bracken_abundance.tsv:** Abundances of Anaerostipus in the samples from bracken
- **inositol_pathway_proteins.faa:** Protein sequencing of the inositol pathway. 
- **inositol_SRR_diamond_0.4_counts.tsv:** Protein hits from diamond analysis.
- **transform.py:** Script to make tables and figures from diamond and bracken data.


# DOWNLOAD data and run DIAMOND
```
#Download all screened fastq files from RAST using a list of RAST ids**
mkdir -p screen.passed.fastq/fasta
cd screen.passed.fastq
cat ../RAST_mgmIDs | parallel -P2 'curl "http://api.metagenomics.anl.gov/1/download/{}?file=299.1" > {}.299.screen.passed.fastq'

#Change file name to their SRR identifiers
sh ../rename_RastID_to_SRR.sh

# Convert fastq to fasta
cat ../SRR_datasets  | parallel -P2 "sed -n '1~4s/^@/>/p;2~4p' {}.299.screen.passed.fastq > fasta/{}.screen.fa"

## Make DIAMOND database from protein file
cd ../
diamond makedb --in inositol_pathway_proteins.faa -d inositol_pathway_proteins

## Run DIAMOND blastx on each read (fasta) file
mkdir diamond_out
cat SRR_datasets | parallel -P1 'diamond blastx -p 1 -e 0.000000001 --db inositol_pathway_proteins -q screen.passed.fastq/fasta/{}.screen.fa -o diamond_out/{}_diamondx_inositol.out'

## Count number of hits in DIAMOND outputs with a cut-off indentity value of at least 40%:
cd diamond_out
awk -F"\t" '{if($3>=40){print $2"\t"FILENAME"\t1" }}' *_diamondx_inositol.out | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | sed 's/_diamondx_inositol.out//g' > ../inositol_SRR_diamond_0.4_counts.tsv
```

# KRAKEN and BRACKEN
***Set up kraken and bracken in your own enviroment and edit locations of files in commands below to make it work***
```
kraken2-build --use-ftp --standard --threads 1 --db KRAKEN2_STANDARD
cat ../SRR_datasets | parallel -P1 'kraken2 --threads 1 --db KRAKEN2_STANDARD screen.passed.fastq/fasta/{}.screen.fa --output {}_kraken2.out --report {}_kraken2.report'
```
**BRACKEN**
```
bracken-build -d ./KRAKEN2_STANDARD/ -t 1 -k 35 -l 350
cat ../SRR| parallel -P1 'python Bracken-2.5/src/est_abundance.py -i {}_kraken2.report -k KRAKEN2_STANDARD/database350mers.kmer_distrib -o {}.report.bracken -t 1'

# Lookup taxonID: 207244 (Genus Anaerostipes) in bracken's est_abundance.py report outputs
grep "207244" *_bracken.report  | awk -F"[_|\t|:]" '{print $1"\t"$(4)"\t"$3}' > Anaerostipes_bracken_abundance.tsv
```

# Generate tables and figures from diamond and bracken output
```
python3 transform.py inositol_SRR_diamond_0.4_counts.tsv Anaerostipes_bracken_abundance.tsv
```
