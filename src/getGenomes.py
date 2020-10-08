#!/usr/bin/env python3

import sys
import argparse
from genetable.entrez import entrez
from genetable.datasets import Datasets
from genetable.utils import dictToPrint


OUTDIR = "genome_out"
SECOL  = ['Organism_Name', 'Assembly_Accession']

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
                        Look and download genomes

    Examples:

        * Look genomes of Yersinia:

            $ getGenomes.py Yersinia
                    
        * Download genomes:

            $ getGenomes.py Yersinia -d

        * Look only the 'Organism_Name' column:

            $ getGenomes.py Yersinia -s Organism_Name

                All available columns are:

                    - Organism_Name
                    - Organism_Kingdom
                    - DefLine
                    - ProjectID
                    - Number_of_Chromosomes
                    - Number_of_Plasmids
                    - Number_of_Organelles
                    - Assembly_Name
                    - Assembly_Accession
                    - AssemblyID
                    - Create_Date
                    - Options

        * Filter 'Organism_Name' names by using list of names:

            $ getGenomes.py Yersinia -c [counter sample file]
""")
    parser.add_argument('term',
                        help='Term')
    parser.add_argument('-c','--counterspps',
                        metavar="",
                        action='store',
                        default=None,
                        help='[Optional] File with counter "Organism_Name"')                        
    parser.add_argument('-d','--download',
                        action="store_true",
                        help='[Optional] if selected, genomes are downloaded')
    parser.add_argument('-s','--columns',
                        action='store',
                        metavar="",
                        nargs= "+",
                        default=SECOL,
                        help='''[Optional] This specify columns. If None, all columns are presented [Default = %s]''' % SECOL)
    parser.add_argument('-o','--outdir',
                        metavar="",
                        action='store',
                        default=OUTDIR,
                        help='[Optional] directory names where genomes will be store [Default = %s]' % OUTDIR)

    args = parser.parse_args()
    return args

def filtercounter(mydict, counterspps):
    
    idx = []
    for n,val in enumerate(mydict['Organism_Name']):
        if not val in counterspps:
            idx.append(n)

    out = {}
    for k,v in mydict.items():
        out[k] = [ v[i] for i in idx ]

    return out

def main():

    # sys.stdout.write("\n")

    opts = getOpts()
    # term = "yersinia"
    if opts.counterspps:
        mycounterlist = []
        with open(opts.counterspps, 'r') as f:
            for i in f.readlines():
                mycounterlist.append(i.strip())

    rawOut = entrez(
                term = opts.term, # opts.term
                db   = "genome" ,
                type = "docsum")

    dictString = rawOut.genomeDS()

    if dictString:
        # rowString  = rawOut._get_type()
        # dictString = genomeDS(rawString)

        if opts.counterspps:
            dictString = filtercounter(dictString, mycounterlist)

        if not dictString:
            sys.stderr.writelines('Empty metadata after filter by countersamples\n')
            sys.stderr.flush()
            exit()

        if opts.download:
            Datasets(dictString, out_dir= opts.outdir).iterate_genome()
            sys.stdout.write("\n")

        else:

            if opts.columns:
                dictString = { c: dictString[c]  for c in opts.columns }

            frows  = dictToPrint(dictString)

            for r in frows:
                print(r)

    else:
        sys.stdout.write("\nNo data for %s\n" % opts.term)
        sys.stdout.flush()    

if __name__ == "__main__":
    main()

