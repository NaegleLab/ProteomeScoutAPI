from proteomeScoutAPI import ProteomeScoutAPI
import proteomeScoutAPI_helpers
import pandas as pd
import argparse



parser = argparse.ArgumentParser(description='Use ProteomeScout API to annotate a dataset that contains peptides and a known accession. Default assumptions are an acc column and a pep column with no modification sites, use optional arguments to set differently.')
parser.add_argument('inputFile', type=str, help='Input File')
parser.add_argument('outputFile', type=str, help='Output File')
parser.add_argument('proteomeScoutFile', type=str, help='Location of ProteomeScout File')
parser.add_argument('-a', '--acc', default='acc', help='Accession column name, if not labeled as acc')
parser.add_argument('-p', '--pep', default='pep', help='Accession column name, if not labeled as acc')
parser.add_argument('-s', '--site', default=False, help='If s or site is TRUE, then find the lowercase PTM in it for centering a peptide and reporting on moifications at a site')
args = parser.parse_args()




inputFile = args.inputFile
outputFile = args.outputFile
proteomeScoutFile = args.proteomeScoutFile
PTM_included = args.site
acc_col = args.acc
pep_col = args.pep

df = pd.read_csv(inputFile)
#First check that the columns that are expected are there
if acc_col not in df:
    raise ValueError("%s not a column for accession found in dataset"%(acc_col))
if pep_col not in df:
    raise ValueError("%s not a column for peptide find in dataset"%(pep_col))

#now instantiate the API
PTM_API = ProteomeScoutAPI(proteomeScoutFile)


for index, row in df.iterrows():
    acc = row[acc_col]
    pep = row[pep_col]

    seq = PTM_API.get_sequence(acc)
    gene_name = PTM_API.get_acc_gene(acc)
    df.at[index, 'gene_name'] = gene_name

    if seq == '-1': 
        #will add -1 to information returned
        continue


    domains = PTM_API.get_domains_harmonized(acc)

    #MODIFICATIONS
    # Return sites and centered peptides for each modification site.
    if PTM_included:
        alignedPeps, seqPosArr, aaArr = proteomeScoutAPI_helpers.returnOrientedPhosphoPeptide(seq, pep)
        pos_aa_arr = []
        for i in range(0, len(seqPosArr)):
            pos_aa_arr.append(aaArr[i]+str(seqPosArr[i]))
        df.at[index, 'modification_sites'] = ';'.join(pos_aa_arr)
        df.at[index, 'aligned_peps'] = ';'.join(alignedPeps)
        #now for each modification, ask if it's in a domain:
        in_domain_arr = []
        for pos in seqPosArr:
            for domain in domains:
                domain_name, domain_start, domain_stop = domain
                if pos >= int(domain_start) and pos <= int(domain_stop):
                    in_domain_arr.append(domain_name)
                    continue #continue back to the next position
        df.at[index, 'site_in_domain'] = ';'.join(in_domain_arr)

        #get Scansite information
        for res_pos in pos_aa_arr:
            scansite = PTM_API.get_Scansite_byPos(acc, res_pos)
            #print("Debug found scanstite")
            #print(scansite)


            if 'scansite_bind' in scansite:
                scansite_bind = scansite['scansite_bind']

                scansite_bind_arr = []
                for bind in scansite_bind:
                    scansite_bind_arr.append(':'.join(bind))
                df.at[index, 'scansite_bind'] = ';'.join(scansite_bind_arr)

            if 'scansite_kinase' in scansite:
                scaniste_kinase = scansite['scansite_kinase']

                scansite_kin_arr = []
                for kin in scaniste_kinase:
                    scansite_kin_arr.append(':'.join(kin))
                df.at[index, 'scansite_kinase'] = ';'.join(scansite_kin_arr)

        #now check if it's been annotated as a known modification site before
        PTMs = PTM_API.get_phosphosites(acc)
        
        found_arr = []
        for pos in seqPosArr:
            found = 0

            for mod in PTMs:
                (pos_mod, res_mod, type_mod) = mod
                #print("Debug, comparing %d to %d"%(pos, int(pos_mod)))
                if pos == int(pos_mod):
                    found = 1
                    #print("FOUND IT!")
                    break

            found_arr.append(str(found))
        df.at[index, 'documented_phosphosite'] = ';'.join(found_arr)





        #print(alignedPeps)
    #find pep in sequence and set as the position of the starting amino acid, then if a PTM, will include all lower cases


    #DOMAINS
    #get domains and the domain strings
    domainStr, archStr = proteomeScoutAPI_helpers.returnDomainArchString(domains)
    df.at[index, 'domains'] = domainStr 
    df.at[index, 'domain_architecture'] = archStr

    #GO TERMS, get them as a string
    try: 
        GO_terms = PTM_API.database[acc]['GO_terms']
        df.at[index, 'GO_terms'] = GO_terms
    except: 
        df.at[index, 'GO_terms'] = -1



df.to_csv(outputFile)
print("Finished.. Wrote new dataframe to outputFile %s"%(outputFile))

