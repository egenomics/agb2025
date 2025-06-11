import pandas as pd

# Load input files 
counts_df = pd.read_csv("tblcounts_asv_melt.csv")  # ASV counts per sample
taxonomy_df = pd.read_csv("tblASVtaxonomy_silva132_v4v5_filter.csv")  # ASV sequences and taxonomy

# Sum ASV counts across all samples 
asv_sums = counts_df.groupby("ASV")["Count"].sum().reset_index()

# Merge with sequences and taxonomy 
merged = asv_sums.merge(taxonomy_df, on="ASV", how="inner")

# Filter out ASVs without genus info 
merged = merged[merged["Genus"] != "<not present>"]

# Keep top 100 ASVs 
top_n = 100  
merged = merged.sort_values(by="Count", ascending=False).head(top_n)

# Normalize counts to total 100,000 reads
merged["NormCount"] = (merged["Count"] / merged["Count"].sum() * 100000).round().astype(int)

# === Write outputs ===

# 1. Write FASTA with only ASV IDs in headers
with open("library_gut.fasta", "w") as fasta:
    for _, row in merged.iterrows():
        fasta.write(f">{row['ASV']}\n{row['Sequence']}\n")

# 2. Write read count table (ASV ID + normalized count)
merged[["ASV", "NormCount"]].to_csv("read_count_gut.tsv", sep="\t", index=False, header=False)

# 3. Save metadata mapping (ASV → Genus, count, sequence)
merged[["ASV", "Genus", "Count", "Sequence"]].to_csv("asv_metadata.tsv", sep="\t", index=False)

print("FASTA written to: library_gut.fasta")
print("Read count file written to: read_count_gut.tsv")
print("Metadata mapping written to: asv_metadata.tsv")
