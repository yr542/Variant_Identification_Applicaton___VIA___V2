import pandas as pd
import numpy as np


# filter the dataFrame (df) by the minimum allele depth (ad) in a particular
# column (name)
def filter_AD(df, name, ad):
    strings=np.array(df[name])
    ADindices=[s.split(":").index("AD") for s in df["FORMAT"]]
    ADs=[strings[i].split(":")[ADindices[i]] for i in range(len(strings))]
    ADs=[int(ad.split(",")[1])/max(int(ad.split(",")[0]),1) for ad in ADs]
    df["AD"]=ADs
    df=df[df["AD"]>ad]
    #print(len(df))
    return df

# filter the dataFrame (df) by minimum depth in a particular column (name)
# if inplace is set to any integer other than 1, it will be filtered into a new data frame
# by default, the function filters in place when the inplace arg is left out of the function call

## Old filter DP removed
def filter_DP(df, name, dp, inplace=1):
    strings = np.array(df[name])
    DPindices = [s.split(":").index("DP") for s in df["FORMAT"]]
    
    DPs = []
    for i in range(len(strings)):
        fixed_str = strings[i].replace(":.", ":0")
        value = fixed_str.split(":")[DPindices[i]]
        DPs.append(int(value))

    df["DP"] = DPs

    if inplace == 1:
        df = df[df["DP"] >= dp].copy()
        del df["DP"]
        return df
    else:
        dfcopy = df[df["DP"] >= dp].copy()
        del dfcopy["DP"]
        return dfcopy


# filter the dataFrame (df) by the maximum number of occurences (cap) of a
# particular zygosity (zyg), e.g. "0/1", in a range of columns
# [namestart,nameend]
def filter_occurences(df, zyg, namestart, nameend, cap):
    freqs=[]
    for i in range(df.columns.get_loc(namestart), df.columns.get_loc(nameend)+1):
        mask=df.iloc[:,i].str.contains(zyg)
        freqs.append(mask)
    freqs=sum(freqs)
    indices=[]
    for key in freqs.keys():
        if freqs[key]>cap:
            indices.append(key)
    df.drop(indices,inplace=True)
    #print(len(df))
    return df

# filter the dataFrame (df) by the maximum population allele frequency (cap)
def filter_AF(df, cap):
    AF_columns=["AF","Kaviar_AF","REGENERON_ALL_AF","gnomad41_genome_AF_grpmax","gnomad41_exome_AF_grpmax"]
    AF_columns = AF_columns + [col + ".1" for col in AF_columns]

    for col in AF_columns:
        if col in df.columns:
            df.loc[df[col]==".",col]="-1"
            df[col]=df[col].astype(float)
            df=df[df[col]<=cap].copy()
    #print(len(df))
    return df

# filter the dataFrame (df) for the zygosity (zyg), e.g. "0/1", in a particular
# column (name)
def filter_zyg(df, name, zyg):
    if name in df.columns:
        df=df[df[name].str.contains(zyg)]
        df=df.drop_duplicates()
    return df

def filter_1x_zyg(df, name, zyg):
    if name in df.columns:
        df=df[df[name].str.startswith(zyg)]
        df=df.drop_duplicates()
    return df

# filter the dataFrame (df) to exclude a certain zygosity (zyg) in a particular
# column (name)
def exclude_zyg(df, name, zyg):
    if name in df.columns:
        df = df[~df[name].str.contains(zyg)]
        df=df.drop_duplicates()
    return df

def exclude_1x_zyg(df, name, zyg):
    if name in df.columns:
        df = df[~df[name].str.startswith(zyg)]
        df=df.drop_duplicates()
    return df
# filter out variants that are "Benign" or "Likely benign"
def filter_benign(df):
    df=df[(df["CLNSIG"].str.contains("enign")==False)]
    return df

# filter the dataFrame (df) for variants with a maximum DP across a list of affected people (names)
# that is greater than the minimum value (dp), a given constant.
# if inplace is 1, it filters df in place; if option is not 1, it filters into a new data frame

# Filter DP Max Changed:

def filter_DP_Max(df, names, dp, inplace=1):
    DPlist = []
    for name in names:
        strings = np.array(df[name])
        DPindices = [s.split(":").index("DP") for s in df["FORMAT"]]
        dp_values = []
        for i in range(len(strings)):
            # Replace any occurrence of ":." with ":0" before splitting
            fixed_str = strings[i].replace(":.", ":0")
            val = fixed_str.split(":")[DPindices[i]]
            dp_values.append(int(val))
        DPlist.append(dp_values)
    DPs = np.max(DPlist, 0)
    df["DP"] = DPs

    if inplace == 1:
        df = df[df["DP"] >= dp].copy()
        del df["DP"]
        return df
    else:
        dfcopy = df[df["DP"] >= dp].copy()
        del dfcopy["DP"]
        return dfcopy


# filter the dataFramd (df) if you only want to keep the rows in which the gene is
# located in a particular chromosome (chrom)
def filter_chr(df, chrom, exclude = False):
    if "Chr" in df.columns:
        if exclude:
            df=df[~df["Chr"].str.contains(chrom)]
        else:
            df=df[df["Chr"].str.contains(chrom)]
    return df

# get which gene in a string of genes (genestring) separated by ;
# is in a list of genes (famgenes), or -1 if none are
def gene_in_list(genestring, famgenes):
    genes = genestring.split(';')
    for i in range(0, len(genes)):
        if genes[i] in famgenes:
            return i
    return -1

# filter the dataframe for only variants in genes associated with the Family
# object's (fam)'s phenotype
def filter_phen(df, fam):
    if len(fam.genes) == 0:
        return pd.DataFrame()

    # get which gene in every gene string is in fam.genes
    gene_locs = [gene_in_list(gene,fam.genes) for gene in df["Gene.refGene"].astype(str)]
    # get a list of booleans:
    # True if a gene in a gene string is in fam.genes, False otherwise
    subset = [num != -1 for num in gene_locs]

    # filter the dataframe by the subset
    df = df[subset]

    # shrink the gene_locs list to only include filtered values
    gene_locs =  [loc for (loc, include) in zip(gene_locs, subset) if include]

    # use the Family's genes-to-n-associated-phenotypes dict to get a list of
    # counts of associated phenotypes for each of the rows.
    # the gene identified by gene_locs is the one passed into the dict.
    counts = [fam.genes[gene.split(";")[loc]] for gene,loc in zip(df["Gene.refGene"],gene_locs)]

    # insert a column containing these counts
    df.insert(3, "phens_matched", counts)

    return df
