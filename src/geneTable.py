#!/usr/bin/env python3

import re
import sys
import argparse
import urllib.request
import matplotlib.pyplot as plt
from collections import Counter
from genetable.entrez import entrez
from genetable.utils import dictToCSV

def getOpts():

    parser = argparse.ArgumentParser(description="Retrieve information from feature table")

    parser.add_argument('string', 
                        metavar='Term',
                        type=str,
                        help='Boolean string which will be used to search sequences')
    parser.add_argument('-q','--qmarkers',
                        metavar="",
                        nargs = "+",
                        default = ["gene","rRNA","tRNA"],
                        help='[Optional] Query markers. String storing markers.[Default =  ["gene" "rRNA" "tRNA"]]')
    parser.add_argument('-p','--plot',
                        action = 'store_true',
                        help='''[Optional] If selected, plot results''')
    parser.add_argument('-y','--hline', metavar="",
                        type = int,
                        default = 3,
                        help='[Optional] Horizontal line in the plot [Default: 3]')
    parser.add_argument('-c','--cutoff', metavar="",
                        default = 10,
                        help='[Optional] Max number of genes at x-axis [Default: 20]')
    parser.add_argument('-r','--cache', metavar="",
                        type=int,
                        default = 100,
                        help='[Optional] Number of sequences downloaded per loop. [Default: 200]')
    parser.add_argument('-o','--out', metavar="",
                        default = None,
                        type= str,
                        help='''[Optional] output file''')
    parser.add_argument('-s','--silent',
                        action = 'store_true',
                        help='''[Optional] If selected, running messages are suppressed''')
    args = parser.parse_args()

    return args

def main():

    args = getOpts()
    # print(args)

    if not args.silent:
        sys.stdout.write("\n")
        sys.stdout.write("\rRetriving information from NCBI...")

    # class calling
    c = entrez(term = args.string,
               type = "ft",
               db   = "nuccore",
               cache = args.cache).feature_table( keyword = args.qmarkers
                                                  ,cutOff = args.cutoff  )
    if c is None:
        print( "Empty feature table under these --qmarkers parameters: %s" %  ", ".join(args.qmarkers) )
        exit()

    if not args.silent:
        sys.stdout.write("\rRetriving information from NCBI...Ok\n")

    outname = '%s_getFeatures' % args.string.replace(" ", "_") if args.out is None else args.out

    if args.plot:

        arr = [i for i in range(0, c.__len__())]

        plt.figure(figsize=(8, 5.5))
        plt.bar(arr
                , c.values()
                , align="center"
                , alpha=0.5)
        plt.xticks(arr
                , c.keys()
                , rotation=87)
        plt.subplots_adjust(bottom=0.33)
        plt.xlabel('Genes')
        plt.ylabel('Frequency')
        plt.title('Gene availability of %s' % args.string)
        plt.axhline(y = args.hline
                    , color = "black")
        plt.savefig(outname + ".png")
        plt.show(block=False)
        plt.close()

        if not args.silent:
            sys.stdout.write("\n")
            sys.stdout.write("Plot stored at %s\n" % (outname + ".png") )

    dictToCSV(dic = c, outname = outname, sep = ",", quiet = args.silent)

if __name__ == "__main__":
    main()