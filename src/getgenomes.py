#!/usr/bin/env python3

import sys
import argparse
# from genetable.entrez import entrez
from genetable.datasets import Datasets

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

                            Get genomes
    Examples:

        * get genomes using a file of accessions:

            $ getgenomes listAcc.txt

            Note: The structure of the `listAcc.txt` file
                  must be like the following:

                  [accession 1],[species name 1]
                  [accession 2],[species name 2]
                  [accession 3],[species name 3]
                   ...         , ...
""")
    parser.add_argument('list',
                        help='list of accesions and species names')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus, mainly used to move files [Default = 1]')                          

    args = parser.parse_args()
    return args

def main():
    opts = getOpts()

    requests = []
    with open(opts.list, 'r') as f:
        for i in f.readlines():
            acc,spps = i.strip().split(",")
            requests.append( (acc, spps) )

    Datasets(
        threads= opts.threads,
        request= requests
    ).iterate_genome()
    
if __name__ == "__main__":
    main()

