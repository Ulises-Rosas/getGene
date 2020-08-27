#!/usr/bin/env python3

import sys
import argparse
from genetable.entrez import entrez
from genetable.datasets import Datasets
from genetable.utils import genomeDS, dictToPrint


OUTDIR = "genome_out"
SECOL  = ['Organism_Name', 'Assembly_Accession']

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='',
                                     epilog="")
    parser.add_argument('term',
                        help='Term')
    parser.add_argument('-l','--lookup',
                        action="store_true",
                        help='[Optional] if selected, only obtain metadata')
    parser.add_argument('-s','--columns',
                        action='store',
                        metavar="",
                        nargs= "+",
                        default=SECOL,
                        help='''if lookup option is used, 
                        this specify columns to be printed.
                        Available columns are: 
                        Organism_Name,
                        Organism_Kingdom, 
                        DefLine, 
                        ProjectID, 
                        Number_of_Chromosomes, 
                        Number_of_Plasmids, 
                        Number_of_Organelles, 
                        Assembly_Name, 
                        Assembly_Accession, 
                        AssemblyID, 
                        Create_Date, 
                        Options
                        [Default = %s]''' % SECOL
                        )
    parser.add_argument('-o','--outdir',
                        metavar="",
                        action='store',
                        default=OUTDIR,
                        help='directory names where genomes will be store [Default = %s]' % OUTDIR 
                        )
    parser.add_argument('-c','--counterspps',
                        metavar="",
                        action='store',
                        default=None,
                        help='Counter sample of species [Default = %s]' % None
                        )
    args = parser.parse_args()
    return args

def main():

    # sys.stdout.write("\n")

    opts = getOpts()
    # term = "yersinia"
    rawOut = entrez(
                term = opts.term, # opts.term
                db   = "genome" ,
                type = "docsum")

    rawString = rawOut._get_type()

    if rawString is not None:
        # rowString  = rawOut._get_type()
        dictString = genomeDS(rawString)

        if opts.lookup:

            if opts.columns:
                dictString = { c: dictString[c]  for c in opts.columns }

            frows  = dictToPrint(dictString)

            for r in frows:
                print(r)

        else:
            Datasets(dictString, out_dir= opts.outdir).iterate_genome()
            sys.stdout.write("\n")

    else:
        sys.stdout.write("\nCheck term: %s\n" % opts.term)
        sys.stdout.flush()    

if __name__ == "__main__":
    main()

