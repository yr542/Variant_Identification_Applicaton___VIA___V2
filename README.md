# Variant Identification Application (VIA) Version 2

Isolates this to a conda environment. VIA automates the process of variant identification, analyzing family data and returning possible candidates based on a range of specific filters and models of inheritance.

## Application Pre-Requisites
- Data Pre-Requisites
  - cleaned data file (output of Mendelian_filtering_WORK.Rmd)
  - pedigree file
  - phenotype file (optional)
- System Pre-Requisites
  - [git](https://git-scm.com/downloads)
  - [python](https://www.python.org/)

Notes on cleaned data file:
- Remove shifted rows from Annovar output (typically 1-10 in large file; check Otherinfo10; should all be .)
- Add headers with sample names (from VCF file)
- Replace all 0|1, 1|1 and 0|0 with 0/1, 1/1 and 0/0


## User Guide

### Installation

To install the application, clone this GitHub repository on your machine. If you have git installed on your computer, you can do this with the following command:

```
git clone https://github.com/yr542/Variant_Identification_Applicaton___VIA___V2.git
```
### Input File Formats

#### Family Pedigree file: csv or txt (tab-delimited)

This file should have the following columns in the following order:

- *Family_ID* (e.g. FIN5)
- *individual_ID* (e.g. FIN5.3)
- *Status* (e.g. Mother, Father, Child, Sibling, Other)
- *Sex* (e.g. Male, Female)
- *Phenotype* (Affected OR Unaffected)

An example would be a tab-delimited table like this:
| Family_ID | individual_ID | Status | Sex    | Phenotype  |
|-----------|---------------|--------|--------|------------|
| FIN1      | FIN1-1        | Father | Male   | Unaffected |
| FIN1      | FIN1-2        | Mother | Female | Unaffected |
| FIN1      | FIN1-3        | Child  | Male   | Affected   |

#### Cleaned Annotated file: csv or txt

This file should be the exact output of the Mendelian_filtering_WORK.Rmd script. Within the defined output format of this script, the important columns used by VIA are:

- *Chr* (Chromosome number)
- *AF_popmax, PopFreqMax, GME_AF, Kaviar_AF, and abraom_freq* (population allele frequencies according to a variety of sources)
- *CLNSIG* (clinical significance of the varient, i.e. benign vs pathogenic)
- *< Individual ID >* (each individual has a column that contains the allelic depth, the zygosity, etc. for each gene)

### Phenotype file: csv or txt (tab-delimited)

This file should have the following columns in the following order:

- *Family_ID* (e.g. FIN5)
- *HPO* (e.g. HP:0000365)

An example would be a tab-delimited table like this:
| Family_ID | HPO                              |
|-----------|----------------------------------|
| FIN1      |HP:0000365                        |
| FIN10     |HP:0001249                        |
| FIN13     |HP:0001249,HP:0012758,HP:0000252  |
| FIN15     |HP:0001251,HP:0001252             |

Notice that when a family has multiple HPO numbers they must be separated by commas (without spaces) in the HPO column.

### Phenotype HPO Ref file: txt (tab-delimited)

This file corresponds to the mapfile input. It should be the latest _phenotype_to_genes.txt_ file from HPO, which can be found at https://hpo.jax.org/app/download/annotation. This file updates every 2 months so the user should check that they have the most updated version of the file in their repository to get the most accurate and current results. Once this file is downloaded it should be in the format compatible with this application without any changes needing to be made.




### Running the Application

The complete application can be run through main.py with the following command:

```
python main.py
```

For greater flexibility, there are also the following optional arguments:
- **_--pedfile_ OR _-p_** : specify the absolute or relative path to the pedigree file. If no argument is specified, the application will look for a file named _Test_Ped.txt_ in the repository's directory
- **_--data_  OR _-d_** : specify the absolute or relative path to the cleaned data file. If no argument is specified, the application will look for a file named _Test_cleaned.txt_ in the repository's directory
- **_--output_ OR _-o_** : specify the file name (including extension and (optionally) the file path) for the ouput file. If no argument is specified, the application will name the output file _filtered.csv_ and place it in the repository's directory
- **_--output_phen_ OR _-op_** : specify the file name (including extension and (optionally) the file path) for the output file when additional phenotype filtering is performed. If no argument is specified, the application will name the output file _filtered_phen.csv_ and place it in the repository's directory
- **_--family_ OR _-f_** : specify a certain family to output an individual csv file for. The default behaviour is to produce a single output file with variants for all families.
- **_--phenfile_ or _-ph_** : specify the absolute or relative path to the phenotype file. If no argument is specified, the application will look for a file named _Test_Phen.txt_ in the repository's directory
- **_--mapfile_ or _-m_** : specify the absolute or relative path to the phenotype-to-gene mapping file. If no argument is specified, the application will look for a file named _phenotype_to_genes.txt_ in the repository's directory. If no such file exists, the user is be prompted to download one.
- **_--nophen_**: specify that no phenotype filtering will be performed.

Any combination of these arguments can be used, and they can be chained together. For example, using all five would look like:

```
python main.py -p <file path> -d <file path> -o <file path> -f <family name> -ph <file path>
```

### Output File Format
  
VIA outputs two csv files with a row for each candidate gene for each individual. In both files the columns are the same as the second input (the cleaned annotated file), except that there are three columns prepended:
  
- *inh model*: The inheritance model(s) (comma-separated) that the variant for that row corresponds to. (e.g. addn; xl,ad; ad, etc.)
- *family*: The Family ID for the individual to whom the variant for that row corresponds to (e.g. FIN5).
- *sample*: The Individual ID for the individual(s) to whom the variant for that row corresponds to (e.g. FIN5.3). Note that the individuals will appear in the same order as their corresponding inheritance models. So if the inheritance models are "xl,ad", and the individuals are "A,B", then the variant is under an xl inheritance model for A and an ad inheritance model for B.'

The first csv file simply has all of the candidate genes for each individual without taking into account the HPO terms. The candidate genes are sorted first in order of sample and then in order of inh model. The second csv file has all of the candidate genes for each individual that also match the phenotype of the affected individuals. The second file also sorts the output in order of family, then in order of the number of HPO terms each candidate 'matches' for a family, then in order of sample. Further, in the second file a column containing the number of HPO terms matched by each candidate is included after the file containing the sample number (so it is the fourth column).

Note: the 'ad, addn' model is not shown for affected individuals with no parents (singletons with or without sibs) in the full output, but will be shown in the 'phen' output.

## Developer Guide

### The Models (models.py)

VIA identifies variants corresponding to four models of inheritance. Note that, for affected singletons, addn and ad are not shown.

#### Model 1: Homozygotes/Hemizygotes (Autosomal Recessive and X-Linked)

Identifies variants in affected individuals who are 1/1 for a given variant. For autosomal recessive, any affected children or siblings are 1/1 and unaffected mothers and fathers are both 0/1. For x-linked, affected siblings and children are 1/1 while unaffected mothers are 0/1 and unaffected fathers are 0/0. The affected person must be male and the variant must be on the x chromosome. For x-linked de novo variants, all of the above hold true except unaffected mothers are 0/0.

#### Model 2: Compound Heterozygotes

Identifies variants in affected individuals who are 0/1 more than or equal to 2 times in a single gene. If the parents are available, variants are only listed when they are correctly phased (one variant comes from one parent and the other from the second parent). If there are more than 2 variants in a single gene in the affected individual, only two of them need to be phased correctly for the variants to be listed.

#### Model 3: De Novo

Identifies variants in affected children/siblings which are not present in either parent. Variants with an allele frequency greater than 0.0005 or a low coverage depth (<6) are removed.

#### Model 4: Autosomal Dominant

Identifies variants in affected children/siblings when _at least one_ of the parents are also affected. Variants with an allele frequency greater than 0.0005 or a low coverage depth (<6) are removed.

### Custom Classes (Family.py)

#### Person

Each sampled individual is a member of the person class. Their characteristic attributes are:
- ID (i.e. FIN5.1)
- sex
- phenotype (unaffected or affected)

#### Family

Each family belongs to the family class. Their characteristic attributes are:
- ID (i.e. FIN5)
- people - a list of Person objects 
- siblings - a list of Person objects whose relationship in the family (relative to the affected individual being studied) is 'sibling'
- mother - a Person object for the mother of the family
- father - a Person object for the father of the family
- child - a Person object for the affected individual being studied

### Custom Filters (filters.py)

Each of these filters are used to pull out candidate variants:
- filter_AD(df, name, ad) - filters a DataFrame (df) by the minimum allele depth (ad) in a paricular column (name)
- filter_DP(df, name, dp, inplace=1) - filters the DataFrame (df) by min depth in a particular column (name). If inplace is set to an integer other than 1, it will filter df into a new data frame, but by default the function filters in place.
- filter_occurences(df, zyg, namestart, nameend, cap) - filters the DataFrame (df) by the max number of occurrences (cap) of a particular zygosity (zyg) in a range of columns
- filter_AF(df, cap) - filters the DataFrame (df) in place by the maximum population allele frequency (cap)
- filter_zyg(df, name, zyg) - filters the DataFrame (df) for the zygosity in a particular column (name)
- exclude_zyg(df, name, zyg) - filters the DataFrame (df) to exclude a certain zygosity (zyg) in a particular column (name)
- filter_benign(df) - filters the DataFrame (df) to exclude variants that are "Benign" or "Likely benign". This filter is not used in any of the models.
- filter_DP_Max(df, names, dp, inplace=1) - filters the DataFrame (df) for variants with a maximum DP across a list of affected people (names) that is greater than the minimum value (dp), a given constant. If inplace is 1, it filters df in place; if it is not, it filters into a new DataFrame
- filter_chr(df, chrom, exclude = False) - filters the DataFrame (df) to keep only the rows in which the gene is located in a particular chromosome (chrom)

## Change Log

### Version 2 (Released June 2025)

This is based on [ischrauwen-lab/variant-filtering](https://github.com/ischrauwen-lab/variant-filtering) which is a fork of [anthonyozerov/variant-filtering](https://github.com/anthonyozerov/variant-filtering). 
* Effectively this is a fork of [ischrauwen-lab/variant-filtering](https://github.com/ischrauwen-lab/variant-filtering), with minor changes. 
* This fork has a docker hub repository available here [yr542/variant_identification_application_via](https://hub.docker.com/r/yr542/variant_identification_application_via).
* The updated stable docker image is placed as a container on Github as well.


