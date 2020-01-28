#!/usr/bin/env python3

import re
import urllib.request
import argparse
import matplotlib.pyplot as plt
from collections import Counter
from entrez import entrez

parser = argparse.ArgumentParser(description="Retrieve information by using API's utils")

parser.add_argument('string', metavar='Term',
                    help='Boolean string which will be used to search sequences')
parser.add_argument('-type', metavar="Retrieval type",
                    default = "fasta",
                    help='<fasta, gb, ...>. Default: fasta')
parser.add_argument('-db', metavar="Database",
                    default = "nuccore",
                    help='<nuccore, sra, ...>. Default: nuccore')
parser.add_argument('--Qmarkers', metavar="Markers",
                    default = "gene,rRNA,tRNA",
                    help='This option only works if Database is "nuccore". String storing markers (e.g. "COI,COX")')
parser.add_argument('-Lmin', metavar="Minimum-Length",
                    default = "",
                    help='This option only works if Database is "nuccore". Minimun length for downloaded sequences')

parser.add_argument('-Lmax', metavar="Maximun-Length",
                    default = "",
                    help='This option only works if Database is "nuccore". Maximun length for downloaded sequences')
parser.add_argument('-ite',
                    action = 'store_true',
                    help='''This option only works if Retrieval type is "gb". This option allows to filter only
                    species-level gg files of both term and group (see below).''')
parser.add_argument('-group', metavar="Iterative",
                    default = "Genus",
                    help='''This option only works if Iterative mode is selected. This option allows to find only
                     species-level "gb" files of a group selected (e.g. Genus). If there was not species-level "gb"
                      files, this value will shift to a higher taxonomic rank (e.g. Family). Default: Genus''')
parser.add_argument('-out', metavar="Iterative",
                    default = "",
                    help='''This option only works if Iterative mode is selected. File name of species-level "gb"
                     files. Default: output.gb''')
parser.add_argument('-ids',
                    action = 'store_true',
                    help='''Get ids of a given request''')
parser.add_argument('--plot',
                    action = 'store_true',
                    help='''plot in some specific functions''')
parser.add_argument('--cutOff', metavar="Markers",
                    default = 20,
                    help='Threshold value')
parser.add_argument('-cache', metavar="Rate",
                    default = 100,
                    help='Number of sequences downloaded per loop. Default: 200')
args = parser.parse_args()

if str(args.db) == "nuccore" and args.type == "ft" and args.plot == True:

    c = entrez(term=str(args.string),
               type=str(args.type),
               db=str(args.db),
               cache = int(args.cache)).feature_table( keyword = str(args.Qmarkers).split(",")
                                               , cutOff= args.cutOff  )

    if c == None:

        print( "Empty feature table under these --Qmarkers parameters: %s" % str(args.Qmarkers) )
    else:

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
        plt.title('Gene availability of %s' % str(args.string))
        plt.axhline(y=3
                    , color="black")
        plt.savefig('%s_GeneAvailability.png' % str(args.string).replace(" ", "_"))
        plt.show(block=False)
        plt.close()

elif str(args.db) == "nuccore" and args.ids == False and args.plot == False:

    entrez(term=str(args.string),
           type=str(args.type),
           db=str(args.db),
           gene_string=str(args.Qmarkers),
           Lmin=str(args.Lmin),
           Lmax=str(args.Lmax),
           printing=True).get_seqs()

elif args.ids:
    for i in entrez(term=str(args.string),
                    db=str(args.db),
                    gene_string=str(args.Qmarkers),
                    Lmin=str(args.Lmin),
                    Lmax=str(args.Lmax))._get_ids():
        print(i)