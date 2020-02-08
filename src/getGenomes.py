#!/usr/bin/env python3

import argparse
from genetable.entrez import entrez
from genetable.utils import genomeDS, dictToPrint

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='',
                                     epilog="")
    parser.add_argument('term',
                        help='Term')
    parser.add_argument('-l','--lookup',
                        action="store_true",
                        help='[Optional] if selected, only obtain metadata')
    parser.add_argument('-o','--out',
                        action='store',
                        default=None,
                        help='Output name. If not stated, results are directly printed at the console [Default = None]'
                        )
    args = parser.parse_args()
    return args

def main():

    opts = getOpts()

    if opts.lookup:
        rawOut = entrez(
                    term = opts.term, # opts.term
                    db   = "genome" ,
                    type = "docsum"
                    )

        if rawOut is not None:
            rowString  = rawOut._get_type()
            dictString = genomeDS(rowString)
            
            frows      = dictToPrint(dictString)
            
            if opts.out is None:
                for r in frows:
                    print(r)
            

if __name__ == "__main__":
    main()

