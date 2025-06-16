import argparse
import pandas as pd
from family import *
from utils import *

if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    argp.add_argument('-p', '--pedfile', default="Test_Ped.txt")
    argp.add_argument('-d', '--data', default="Test_cleaned.txt")
    argp.add_argument('-o', '--output', default="filtered.csv")
    argp.add_argument('-op', '--output_phen', default="filtered_phen.csv")
    argp.add_argument('-f', '--family', default="")
    argp.add_argument('-ph', '--phenfile', default="Test_Phen.txt")
    argp.add_argument('-m', '--mapfile', default="phenotype_to_genes.txt")
    argp.add_argument('--nophen', default = False, action = 'store_true')

    args = argp.parse_args()

    # get a dict of families from the pedfile
    families = get_families(args.pedfile)

    
    if not args.nophen:
        print("Getting relevant genes for family phenotypes...")
        # give each family a list of genes relevant to their phenotype
        load_phen(families, args.phenfile, args.mapfile)

    # read in the file containing variants
    df = pd.read_csv(args.data, sep='\t', low_memory=False)
    #check that there are no errors, and remove rows with errors.
    df = verify(df)

    # csv with variants in one family
    if args.family != "":
        fam = families[args.family]
        fam_variants = df.copy()
        for person in fam.people:
            filt = filter_zyg if person.phen == "Unaffected" else exclude_zyg
            fam_variants = filt(fam_variants, person.ID, "0/0")
        fam_variants.to_csv(fam.ID + ".csv")


    # empty dataframes for results with and without phenotype filter
    result = pd.DataFrame()
    result_p = pd.DataFrame()

    for fam in families.values():

        print("Filtering", fam.ID + '...')

        # get a dataframe of variants for the family,
        # without phenotype filter
        famresult = filter_family(df, fam, phenfilter = False)
        # append it to the results
        result = pd.concat([result,famresult])

        if not args.nophen:
            # get a dataframe of variants for the family,
            # with phenotype filter
            famresult_p = filter_family(df, fam, phenfilter = True)
            # append it to the results
            result_p = pd.concat([result_p,famresult_p])

    # organize result first by sample and then by inh model
    result = result.sort_values(['sample', 'inh model'])
    #save result
    result.to_csv(args.output)
    print(result)

    #save result with phenotype filter
    if not args.nophen:
        result_p = result_p.sort_values(['family','phens_matched','sample'])
        result_p.to_csv(args.output_phen)
        print(result_p)
