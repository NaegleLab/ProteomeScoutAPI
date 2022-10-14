

# Common Helper functions to move in and out of strings/tuples from ProteomeScout, align peptides around modifications, etc.


def returnDomainArchString(domains):
    """
    Given a list of domains, create a printable string for them

    Parameters
    ----------
    domainList: list of tuples
        What is returned from PTM_API.get_domains*

    Returns
    -------
    domainString: string
        A string with domains separated by ';' and as domainName:start:stop
    domainArchString: string
        A string with the architecture, ordered by how they appear in protein and separated by '-', e.g. domain1-domain2-domain3

    """
    domainStruct = {}
    domainList = []
    for domain_entry in domains:
        name, start, stop = domain_entry
        #print("DEBUG: found %s"%(name))
        domainList.append("%s:%s:%s"%(name, start, stop))
        domainStruct[int(start)] = name
    domainStarts = list(domainStruct.keys())

    archList = []
    for domain_start in domainStarts:
        archList.append(domainStruct[domain_start])


    domainStr = ';'.join(domainList)
    archStr = '~'.join(archList)

    return domainStr, archStr


def returnGOStrings(GO_terms):
    """
    Given a GO terms from the API, return the strings for each of the key types of GO terms.
    
    """


def returnStartOfPeptidePosition(proteinSeq, pep):
    """
    Return the protein position that is the start of the peptide. Will return -1 if peptide not found
    """
    pos_start = -1

    try: 
        indStart = proteinSeq.index(pep.upper())+1
    except ValueError:
        print("ERROR: pep %s not found"%(pep))
    return pos_start

def returnOrientedPhosphoPeptide(proteinSeq, pep, flank=7):
    """
    For every lowercase amino acid in a pep sequence, match it to 
    protein and return an oriented peptide, with one lc letter in center,
    and flanked by x-amino acids. This currently looks only for s, t, and y. Can update find_phospho to find other modification types

    Parameters
    ----------
    proteinSeq: str
        PTM_API.get_sequence(acc) - protein sequence
    pep: str
        peptide sequence to match to protein sequence
    flank: int
        default 7. The number of amino acids to add to the n and c terminal side of the phosphoryaltion site.

    """
    pTyr_posArr = find_phospho(pep)
    seqPosArr = []
    alignedArr = []
    aaArr = []
    try: 
        indStart = proteinSeq.index(pep.upper())
    except ValueError:
        print("ERROR: pep %s not found"%(pep))
    else:
        #
        for pTyr_pos in pTyr_posArr:
            seq_pos = indStart+pTyr_pos
            seqPosArr.append(seq_pos+1)
            aligned = proteinSeq[seq_pos-flank:seq_pos]+proteinSeq[seq_pos].lower()+proteinSeq[seq_pos+1:seq_pos+flank+1]
            #wrote this assuming we weren't at C/Nterminal for any of these
            alignedArr.append(aligned)
            aaArr.append(proteinSeq[seq_pos])
            #check if it is know to be phosphorylated
            
    return alignedArr, seqPosArr, aaArr

def find_phospho(s):
    return [i for i, letter in enumerate(s) if (letter == 't' or letter=='s' or letter=='y')]


