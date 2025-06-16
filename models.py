# This file is for functions that apply the filters to the families based on
# the different inheritance models (ar, xl, xldn, ch, addn, ad)
# where:
# ar: autosomal recessive (model #1)
# xl: x-linked (model #1)
# xldn: x-linked de novo (model #1)
# ch: compound heterozygous (model #2)
# addn: autosomal dominant de novo (model #3)
# ad: autosomal dominant (model #4)

# The functions should take in the variant data frame and a Family object,
# and output a data frame of possible variants.

from family import Family
from filters import *
import pandas as pd

#add_columns adds three columns to the dataFrame df containing info to be outputted for
#candidate variants for a specific Family object fam and the relevant inheritance
#model number modelno
def add_columns(df, fam, model):
    df.insert(0, "inh model", model)
    df.insert(1, "family", fam.ID)
    df.insert(2, "sample", fam.child.ID)

# ad_model takes in a data frame and Family object and returns a new data frame
# containing candidate variants
def ad_model(df, fam, include_singleton = False):
    min_allelic_depth = 6  # will filter for 6x coverage minimum for at least one affected individ
    numAffected = 0
    newdf = df.copy()
    newdf = filter_AF(newdf, .0005)  # filters all AF cols for entries <= .0005
    dpdf = pd.DataFrame()
    names = []
    for person in fam.people:
        if person.affected:
            names.append(person.ID)
            numAffected += 1
            newdf = filter_zyg(newdf, person.ID, "0/1")  # filters for 0/1 entries for affected individs
        else:
            newdf = filter_zyg(newdf, person.ID, "0/0")  # filters for 0/0 entries for unaffected individs

    # returns an empty Data Frame if nothing should be output for this model (<= 1 affected individs
    # or they are a singleton)
    noparents = not fam.hasFather and not fam.hasMother

    if numAffected == 0: return pd.DataFrame()
    elif not include_singleton and noparents: return pd.DataFrame()
    else:
        newdf = filter_DP_Max(newdf, names, min_allelic_depth,0)
        add_columns(newdf, fam, "ad")  # adds on columns with family info
        return newdf

# de_novo_model takes a dataframe (the cleaned data) and a family object
# return value: a new dataframe with all possible de novo candidate
# genes
def de_novo_model(df, fam, include_singleton = False):

    # re-filter for MAF
    revised_df = df.copy()
    revised_df = filter_AF(revised_df, .0005)
    
    # keep track of number of individuals we are identifying variants
    # for
    num_affected = 0

    # If either mother or father is affected, no de novo, so return
    # empty data frame
    if fam.mother.affected or fam.father.affected:
        return pd.DataFrame()

    noparents = not fam.hasMother and not fam.hasFather
    if noparents and not include_singleton:
        return pd.DataFrame()
   
    # filter child for all 0/1
    if fam.child.ID != "":
        num_affected += 1
        revised_df = filter_zyg(revised_df, fam.child.ID, "0/1")
        
	# filter to make sure DP is at least 6x
        revised_df = filter_DP(revised_df, fam.child.ID, 6)

    # filter parents for 0/0
    revised_df = filter_zyg(revised_df, fam.father.ID, "0/0")
    revised_df = filter_zyg(revised_df, fam.mother.ID, "0/0")

    # filter siblings to identify more candidate genes
    for sib in fam.siblings:
        if sib.affected:
            num_affected += 1
            revised_df = filter_zyg(revised_df, sib.ID, "0/1")
            revised_df = filter_DP(revised_df, sib.ID, 6)
        else:
            revised_df = filter_zyg(revised_df, sib.ID, "0/0")
    
    if num_affected:

        # add on the columns with family info
        add_columns(revised_df, fam, "addn")
        return revised_df

    # if no affected individuals, return empty data frame
    return pd.DataFrame()

def cmpd_het_model(df, fam):
   
    # keep track of individuals we are identifying variants for
    num_affected = 0

    # filter child for having 0/1 in >=2 variants of the same gene
	
    # create newdf to include all instances of child 0/1
    if fam.child.ID != "":
        num_affected += 1
        newdf = df.copy()
        newdf = filter_zyg(newdf, fam.child.ID, "0/1")
	
        # use Gene.refGene column to create new column "Gene' with
        # gene names (deals with semicolon issue in some genes)
        partitioned_gene=0
        newdf.loc[:,'Gene'] = ''
        newdf=newdf.dropna(subset=["Gene.refGene"])
        newdf = newdf.reset_index(drop=True)
        newdf['Gene'] = newdf["Gene.refGene"].copy().str.partition(";")[0]
        
        # create new df where genes must meet the following criteria: 
	# there must be at least 1 0/1 variant in mother that is not in father
	# and at least 1 0/1 variant in father that is not in mother
        finaldf = newdf[newdf.duplicated(subset=['Gene'], keep=False)] 
        if(fam.hasFather or fam.hasMother):
            genes = finaldf["Gene"].unique()
            both_available = fam.hasFather and fam.hasMother
            if not both_available:
                available_ID = fam.father.ID if fam.hasFather else fam.mother.ID
            for gene in genes:
                genedf = finaldf[finaldf["Gene"]==gene]
                if(both_available):
                    mom = sum(genedf[fam.mother.ID].str.contains("0/1") &
                              genedf[fam.father.ID].str.contains("0/0"))
                    dad = sum(genedf[fam.mother.ID].str.contains("0/0") &
                              genedf[fam.father.ID].str.contains("0/1"))
                    if(mom==0 or dad==0):
                        finaldf = finaldf[finaldf["Gene"]!=gene]
                else:
                    parentvariants = sum(genedf[available_ID].str.contains("0/1"))
                    parentnonvariants = sum(genedf[available_ID].str.contains("0/0"))
                    if(parentvariants == 0 or parentnonvariants == 0):
                        finaldf = finaldf[finaldf["Gene"]!=gene]

        # delete the gene column we created
        del finaldf['Gene']
        
    # add on the columns with family info
    if num_affected:
        add_columns(finaldf, fam, "ch")
        return finaldf

def xl_model(df, fam):
    newdf = df.copy()
    x_df = filter_chr(newdf, "chrX")
    for person in fam.people:
        if person.affected:
            if not person.male:
                return pd.DataFrame()
            x_df_1 = filter_zyg(x_df, person.ID, "1/1")
            x_df_2 = filter_1x_zyg(x_df, person.ID, "1:")
            x_df = x_df_1.append(x_df_2)
        if person.unaffected:
            x_df_1 = exclude_zyg(x_df, person.ID, "1/1")
            x_df_2 = exclude_1x_zyg(x_df, person.ID, "1:")
            x_df = x_df_1.append(x_df_2)
            if person.male:
                x_df_1 = filter_zyg(x_df, person.ID, "0/0")
                x_df_2 = filter_1x_zyg(x_df, person.ID, "0:")
                x_df = x_df_1.append(x_df_2)
    if not fam.father.affected:
        x_df = filter_zyg(x_df, fam.mother.ID, "0/1")

    add_columns(x_df, fam, "xl")
    return(x_df)

def xldn_model(df, fam):
    newdf = df.copy()
    x_df = filter_chr(newdf, "chrX")
    if fam.mother.affected or fam.father.affected:
        return pd.DataFrame()
    if fam.child.female:
        return pd.DataFrame()
    # filter child for all 0/1
    for person in fam.people:
        if person.affected:
            if not person.male:
                return pd.DataFrame()
            x_df_1 = filter_zyg(x_df, person.ID, "1/1")
            x_df_2 = filter_1x_zyg(x_df, person.ID, "1:")
            x_df = x_df_1.append(x_df_2)
        if person.unaffected:
            x_df_1 = filter_zyg(x_df, person.ID, "0/0")
            x_df_2 = filter_1x_zyg(x_df, person.ID, "0:")
            x_df = x_df_1.append(x_df_2)
    add_columns(x_df, fam, "xldn")
    return(x_df)
    

# ar_model takes the data frame and Family object. Returns: a new data frame containing
# all possible autosomal recessive candidate genes
def ar_model(df, fam):
    min_allelic_depth = 6  # will filter for 6x coverage minimum for at least one affected individ
    numAffected = 0
    newdf = df.copy()
    newdf = filter_AF(newdf, .005) 
    newdf = filter_chr(newdf, "chrX", exclude=True)
    dpdf = pd.DataFrame()
    names = []
    for person in fam.people:
        if person.affected:
            names.append(person.ID)
            numAffected += 1
            newdf = filter_zyg(newdf, person.ID, "1/1") 
        else:
            newdf = exclude_zyg(newdf, person.ID, "1/1")
            if person == fam.father or person == fam.mother:
                newdf = filter_zyg(newdf, person.ID, "0/1")
    # returns an empty Data Frame if nothing should be output for this model (<= 1 affected individs)
    if numAffected < 1:
        return pd.DataFrame()
    else:
        newdf = filter_DP_Max(newdf, names, min_allelic_depth,0)
        add_columns(newdf, fam, "ar")  # adds on columns with family info
        return newdf
    
