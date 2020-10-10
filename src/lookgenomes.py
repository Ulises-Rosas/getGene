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

                            Look at genomes metadata

    Examples:

        * Look genomes of Yersinia:

            $ lookgenomes Yersinia
                    
        * Filter 'Organism_Name' names by using list of names:

            $ lookgenomes Yersinia -c [counter sample file]

            Note: Counter samples is a list of non-target 
                  species in plain text
""")
    parser.add_argument('term',
                        help='Term')
    parser.add_argument('-c','--counterspps',
                        metavar="",
                        action='store',
                        default=None,
                        help='[Optional] File with counter samples')
    args = parser.parse_args()
    return args

def printthisout(results):

    if results:
        sys.stdout.write(
            "\t".join([
                    'accession'     ,
                    'seq_length'    ,
                    'conting_n50'   ,
                    'sci_name'      ,
                    'assembly_level',
                    'title'
                ]) + "\n"
            )
        sys.stdout.flush()

        for tmpdict in results:
            sys.stdout.write(
                "\t".join([
                        tmpdict['accession']       ,
                        tmpdict['seq_length']      ,
                        str(tmpdict['conting_n50']),
                        tmpdict['sci_name']        ,
                        tmpdict['ass_level']       ,
                        tmpdict['title']
                    ]) + "\n"
                )
            sys.stdout.flush()


def filtercounter(rawOut, counterspps):
    # counterspps = ['Denticeps clupeoides']
    out = []

    for mydict in rawOut:
        tmpspps = mydict['sci_name']

        if not tmpspps in counterspps:
            out.append(mydict)
            
    return out

def main():
    opts = getOpts()
    mycounterlist = []

    if opts.counterspps:
        with open(opts.counterspps, 'r') as f:
            for i in f.readlines():
                mycounterlist.append(i.strip())

    # datsetsclass = Datasets()
    rawOut = Datasets().taxon_descriptor(spps=opts.term)

    if not rawOut:
        exit()

    if mycounterlist:
        rawOut = filtercounter(rawOut, mycounterlist)

    printthisout(rawOut)

    
if __name__ == "__main__":
    main()

