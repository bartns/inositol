#/usr/bin/python3
import sys
import numpy as np 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

diamond_hit_counts = sys.argv[1]
bracken_data_file = sys.argv[2]

diamond_counts = open(diamond_hit_counts,"r").readlines()
# create Pandas DataFrame from diamond hit counts
df = pd.DataFrame()
for i,line in enumerate(diamond_counts):
    protein = line.split("\t")[0]
    sample = line.split("\t")[1]
    count = int(line.strip().split("\t")[2])
    if not protein in df.index.values: 
        df = df.append(pd.Series(name=protein,dtype=int))
    if not sample in df: 
        df[sample] = 0
    df[sample][protein] = count
# replace NaN values with 0
df = df.fillna(0)
# sort columns 
df = df.reindex(sorted(df.columns), axis=1)
df.to_csv(r"pathway_counts_table.csv")
#print(df, sep='\n')

# calculate normalized pathway abundance per sample. 
# method: Vital et al. DOI: 10.1128/mBio.00889-14
totalcounts = {}
for line in open("totalReads_samples.tsv","r").readlines(): 
    totalcounts[line.split()[0]] = int(line.strip().split()[1])

abundance_df = pd.DataFrame(columns=["#",'sample','abundance'],index=["#"])
for i,sample in enumerate(df):
    # 30249 = total basepairs in inositol pathway
    # grep -v ">" inositol_pathway_proteins.faa | wc | awk '{print ($3-$1)*3}'
    sample_abundance = df[sample].sum()/(totalcounts[sample]*30249/4000000)
    abundance_df = abundance_df.append({'#':i,"sample":sample,"abundance":sample_abundance}, ignore_index=True)

abundance_df.set_index('#', inplace = True)
abundance_df.sort_values(by=['sample'])
abundance_df.to_csv("pathway_abundance.csv")

# make DataFrame binary (1 when cell > 0)
bin_df = df.copy()
for sample in bin_df: 
    bin_df[sample] = np.asarray(np.ceil(bin_df[sample] / df[sample].sum()), dtype=int)
bin_df.to_csv(r"pathway_absence-presence_table.csv")

# Absence presence heatmap (from binary matrix)
sns.set()
plt.figure(figsize=(35, 8))
heatmap = sns.heatmap(bin_df, linewidths=.5,annot=False,cmap=sns.light_palette("steel blue", input="xkcd"))
fig = heatmap.get_figure()
fig.savefig("pathway_protein_presence-absence_heatmap.svg",format='svg', dpi=300)

# Pathway Abundance barplot
plt.figure(figsize=(24, 2))
barplot = sns.barplot(x='sample',y='abundance',data=abundance_df, color="steelblue")
axes = barplot.axes
axes.set_ylim(0,1)
axes.set_yticks([0,0.25,0.50,0.75,1.00])
fig = barplot.get_figure()
fig.savefig("pathway_normalized-abundance_barplot.svg",format='svg', dpi=300)

## Create Bracken abundance DataFrame
bracken_data_file = sys.argv[2]
bracken_abundance = open(bracken_data_file,"r").readlines()
bracken_df = pd.DataFrame(columns=["#",'sample','abundance'],index=["#"])
for i,line in enumerate(bracken_abundance):
    sline = line.strip().split("\t")
    bracken_df = bracken_df.append({'#':i,"sample":sline[0],"abundance":sline[2]}, ignore_index=True)
bracken_df.set_index('#', inplace = True)
bracken_df.sort_values(by=['sample'])

print(bracken_df, sep='\n')

# Bracken abundance in Barplot
plt.figure(figsize=(24, 2))
barplot = sns.barplot(x='sample',y='abundance',data=bracken_df, color="steelblue")
axes = barplot.axes
axes.set_ylim(0,1)
axes.set_yticks([0,0.25,0.50,0.75])
fig = barplot.get_figure()
fig.savefig("pathway_normalized-abundance_barplot.svg",format='svg', dpi=300)


    
