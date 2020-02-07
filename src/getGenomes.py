#!/usr/bin/env python3

import argparse
from genetable.entrez import entrez

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='',
                                     epilog="")
    parser.add_argument('term',
                        help='Term')
    parser.add_argument('-l','--look',
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

    if opts.look:
        rawOut = entrez(
                    term = opts.term,
                    db   = "genome" ,
                    type = "docsum"
                    )

        if rawOut is not None:
            print( rawOut._get_type() )        

if __name__ == "__main__":
    main()

