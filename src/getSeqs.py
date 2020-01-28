#!/usr/bin/env python3
import argparse
from  entrez import entrez

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='',
                                     epilog="")
    parser.add_argument('term',
                        metavar='term',
                        default=None,
                        help='Term')
    parser.add_argument('--db',
                        metavar="str",
                        default="nuccore",
                        help='[Optional] entrez database to look for [Default = "nuccore"]')
    parser.add_argument('--rettype',
                        metavar="str",
                        default="fasta",
                        help='[Optional] retrieve type [Default = "fasta"]')
    parser.add_argument('--cache',
                        metavar="int",
                        default=200,
                        help='[Optional] numbers of sequences to download per iteration [Default = 200]')
    parser.add_argument('--Lmin', '-l',
                        metavar ="min length",
                        default = "",
                        help = 'This option only works if database is "nuccore". Minimun length for downloaded sequences')
    parser.add_argument('--Lmax', '-L',
                        metavar = "max length",
                        default = "",
                        help = 'This option only works if Database is "nuccore". Maximun length for downloaded sequences')
    parser.add_argument('--idsOnly', '-d',
                        action="store_true",
                        help = 'If selected, only IDS of the request will be the outcome')
    parser.add_argument('--out', '-o', metavar="str",
                        action='store',
                        default=None,
                        help='Output name. If not stated, results are directly printed at the console [Default = None]'
                        )
    args = parser.parse_args()
    return args

def printOutput(out, file):

    if file is None:
        for i in out:
            print(i)
    else:
        f = open(file, "w")
        for i in out:
            f.write(i + "\n")
        f.close()

def main():
    opts = vars(getOpts())

    # opts = {
    #     "Lmax":'',
    #     "Lmin":'',
    #     "cache":200,
    #     "db":'nuccore',
    #     "idsOnly":True,
    #     "out":None,
    #     "rettype":'fasta',
    #     "term":'JL597290:JL597291[ACCN]'
    # }
    # print(opts)

    classEntry = entrez(term  = opts["term"],
                        type  = opts["rettype"],
                        db    = opts["db"],
                        Lmin  = opts["Lmin"],
                        Lmax  = opts["Lmax"],
                        cache = opts["cache"],
                        printing = True if opts["out"] is None else False
                        )

    if opts["idsOnly"]:
        out = classEntry._get_ids()
        printOutput(out, opts["out"])
        exit()

    if opts["out"] is None:
        classEntry.get_seqs()
    else:
        printOutput(out, opts["out"])

if __name__ == "__main__":
    main()