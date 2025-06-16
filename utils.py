import pandas as pd
from family import *
from models import *

# get a dict of family IDs as keys and Family objects as values
# from the PED file (pedfile)
def get_families(pedfile):

    # read the pedfile into a dataframe
    pedDf = pd.read_csv(pedfile, sep='\t')

    #create an empty dict to store the families
    families = {}
    for i in range(0, len(pedDf)):

        Fam_ID = pedDf["Family_ID"][i]
        # if we have not encountered this family yet,
        # create a new Family object, else get the already-created
        # one from the dict
        if Fam_ID not in families.keys():
            fam = Family(Fam_ID)
            families[Fam_ID] = fam
        else:
            fam = families.get(Fam_ID)

        # create a Person object representing this row of the
        # PED dataframe
        ID = pedDf["individual_ID"][i]
        status = pedDf["Status"][i]
        sex = pedDf["Sex"][i]
        phen = pedDf["Phenotype"][i]
        newperson = Person(ID, sex, phen)
        # add the Person to the Family's list of people
        fam.people.append(newperson)

        # update variables in the Family object based on this
        # Person's position in the famile
        if status == "Father":
            fam.father = newperson
            fam.hasFather = True
        elif status == "Mother":
            fam.mother = newperson
            fam.hasMother = True
        elif status == "Child":
            fam.child = newperson
        elif status == "Sibling":
            fam.siblings.append(newperson)

    return families

import os
# give each family in the list of families (families) a list of genes
# relevant to their phenotype.
# the phenotype is taken from the phenotype file (phenfile) and the mapping
# from HPO number to genes is taken from (mapfile) or, if it does not exist,
# is downloaded.
def load_phen(families, phenfile, mapfile):

    # if the mapfile does not exist in the current directory
    if not os.path.isfile(mapfile):

        # offer to download it
        answer = input("No phenotype-to-gene mapping found. Download one? [y/N]: ").lower()

        if answer=="y":
            # download the mapfile
            import urllib.request
            url = "https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt"
            print("Downloading now. It might take a while...")
            urllib.request.urlretrieve(url,mapfile)
            print("Finished downloading.")
        else:
            print("exiting now")
            exit()

    # read the mapfile into a dataframe
    phen_to_genes = pd.read_csv(mapfile, sep = '\t', header = None, comment = '#')

    # rename the columns of the dataframe
    phen_to_genes.columns = ["HPO-id", "HPO label", "gene-id",
            "gene-symbol", "additional-info",
            "source", "disease-ID"]

    # read the phenfile into a dataframe
    phenDf = pd.read_csv(phenfile, sep='\t')
    for i in range(0, len(phenDf)):

        # get the family ID in this row of phenotype dataframe
        Fam_ID = phenDf["Family_ID"][i]
        if Fam_ID in families:
            # get the Family object that corresponds to this ID
            fam = families.get(Fam_ID)
            # get the HPO numbers from this row of the phenotype dataframe
            hpo = phenDf["HPO"][i]
            # get a list of HPO numbers by splitting them by the comma
            fam.HPO = hpo.split(',')
            # for each HPO number,
            for HPO in fam.HPO:
                # use the phenotype-to-gene dataframe to get a list of genes
                # associated with that HPO number
                genes = phen_to_genes[phen_to_genes["HPO-id"]==HPO]['gene-symbol'].tolist()
                # get a list of the number of phenotypes associated with each gene
                # (it will be just 1 if we have not encountered this gene
                #  in this family yet, and the existing number + 1 otherwise)
                gene_nums = [fam.genes[gene]+1 if gene in fam.genes else 1 for gene in genes]
                # update the family's dict of associated genes and their
                # numbers of phenotypes
                fam.genes.update(dict(zip(genes, gene_nums)))

# Checks that DP is in every row in the FORMAT column
def verify(df):
    # get boolean series, with True if a row is bad
    badrows = ~df["FORMAT"].str.contains("DP")

    if any(badrows):
        print("The following rows do not contain 'DP' in their 'FORMAT' column:")
        # print indices as they would appear in spreadsheet software
        print([idx + 2 for idx in df[badrows].index])
        # remove the bad rows from the df
        df = df[~badrows]
        print("They will not be considered.")

    return df

# generate a list of subfamilies centered on every affected individual in the
# Family object (fam)
def generate_subfamilies(fam):

    # the original Family is a subfamily
    subfamilies = [fam]

    for person in fam.people:
        if person.affected and person != fam.child:
            # the family ID is the same
            subfamily = Family(fam.ID)

            # the child of this family is the affected individual
            subfamily.child = person

            # if the person is a parent in the original family, the
            # opposite person is the other parent
            opposite = ""
            if person == fam.father or person == fam.mother:
                opposite = fam.mother if person == fam.father else fam.father

            # if the person is a sibling in the original family,
            # then we give this subfamily the same parents,
            # and a list of siblings including all of the siblings
            # in the original family that are not this person, as well
            # as the child in the original family
            if person in fam.siblings:
                if fam.hasFather:
                    subfamily.father = fam.father
                    subfamily.hasFather = True
                if fam.hasMother:
                    subfamily.mother = fam.mother
                    subfamily.hasMother = True
                newsibs = fam.siblings.copy()
                newsibs.remove(person)
                subfamily.siblings = newsibs + [fam.child]

            # add all the people in the original family who are not the
            # opposite parent to the list of people in the subfamily
            for p in fam.people:
                if p != opposite:
                    subfamily.people.append(p)
            subfamilies.append(subfamily)
    return subfamilies

# combine multiple instances of the same variant in df into one row,
# getting a new dataframe
def combine_duplicates(df):
    # use a generated location string to determine if two variants
    # are the same or not. These location strings should be different
    # if and only if the two variants are different
    df["loc"] = [chrom + str(start) + str(end) for chrom, start, end in
                        zip(df['Chr'], df['Start'], df["End"])]

    # get a list of the unique location strings
    uniquelocs = df["loc"].unique()

    # create an empty dataframe to store the output
    combined = pd.DataFrame()

    # for every unique location string,
    for loc in uniquelocs:
        # get a dataframe of the rows in df with that location string
        rows = df[df["loc"] == loc]
        # since the rows are identical except for the sample and inh model
        # in formation, get a dataframe containing just the first of the rows
        output = rows.head(1).copy()

        # create a string that is all of the samples in the rows separated
        # by commas
        samplestring = ""
        for sample in rows["sample"]:
            samplestring += sample + ","
        samplestring = samplestring[:-1]  # remove last comma

        # create a string that is all of the inheritance models in the rows
        # separated by commas
        modelstring = ""
        for model in rows["inh model"]:
            modelstring += model + ","
        modelstring = modelstring[:-1]  # remove last comma

        # put these strings into the output row and delete the column
        # containing the location string
        output["sample"] = [samplestring]
        output["inh model"] = [modelstring]
        del (output["loc"])

        #add the output row to the dataframe
        combined = pd.concat([combined, output])

    return combined

# filter a dataframe of variants (df), getting the ones for which
# inheritance models for the Family object (fam) fit.
# apply the phenotype filter if phenfilter is True
def filter_family(df, fam, phenfilter):

    # generate a list of subfamilies centered on each affected individual
    subfamilies = generate_subfamilies(fam)

    # empty dataframe for results in this family
    famresult = pd.DataFrame()

    if phenfilter and len(fam.genes) == 0:
        print("Warning: no phenotypes listed for", fam.ID, "in the phenotype file.")
        return famresult

    # add model results for each subfamily
    for subfam in subfamilies:
        famresult = pd.concat([famresult, ad_model(df, subfam, include_singleton = phenfilter)])
        famresult = pd.concat([famresult, ar_model(df, subfam)])
        famresult = pd.concat([famresult, xl_model(df, subfam)])
        famresult = pd.concat([famresult, xldn_model(df, subfam)])
        famresult = pd.concat([famresult, de_novo_model(df, subfam, include_singleton = phenfilter)])
        famresult = pd.concat([famresult, cmpd_het_model(df, subfam)])
    # combine multiple instances of the same variant into one row
    famresult = combine_duplicates(famresult)

    # additionally apply the phenotype filter if requested
    if phenfilter:
        famresult = filter_phen(famresult, fam)

    # print helpful output
    phenfilterstring = 'with   ' if phenfilter else 'without'
    number = '{0:4d}'.format(len(famresult))
    print(number,'candidates',phenfilterstring,'phenotype filter')

    # return the result
    return famresult
