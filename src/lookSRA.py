import sys
import re
import argparse
from genetable.entrez import entrez
# from genetable.datasets import Datasets
# from genetable.utils import genomeDS, dictToPrint
import xml.etree.ElementTree as ET


def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
                        Look SRA information
    Examples:

        * Get metadata from a given term:

            $ lookSRA.py 'Yersinia intermedia'

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

def formatmytype(mytype, countersamples):

    tree = ET.fromstring(mytype)

    out= []

    for child in tree:

        expxml = child[1].text
        runa   = child[2].text

        if expxml:

            libsrcpatt = '.*<LIBRARY_SOURCE>(.+)</LIBRARY_SOURCE>.*'
            srcname  = re.sub(libsrcpatt, "\\1", expxml)

            libstrpatt = '.*<LIBRARY_STRATEGY>(.+)</LIBRARY_STRATEGY>.*'
            straname = re.sub(libstrpatt, "\\1", expxml)

            platpatt = '.*instrument_model="(.*?)".*'
            platname = re.sub(platpatt, "\\1", expxml)

            scipatt = '.*ScientificName="(.*?)".*'
            sciname = re.sub(scipatt, "\\1", expxml)

            if countersamples:
                if sciname in countersamples:
                    expxml = None

        if runa:
            accpatt = '.* acc="(.*?)".*'
            runacc  = re.sub(accpatt, "\\1", runa)

            ispupatt = '.* is_public="(.*?)".*'
            ispublic = re.sub(ispupatt, "\\1", runa)
        
        if expxml and runa:
            tmp = {
                'sciname'   : sciname,
                'platname'  : platname,
                'runacc'    : runacc,
                'ispublic'  : ispublic,
                'strategy'  : straname,
                'lib_source': srcname
            }
            out.append(tmp)

    return out

def printthisout(results):

    if results:
        sys.stdout.write(
            "\t".join([
                    'accession'   ,
                    'species_name',
                    'platform'    ,
                    'is_public'   ,
                    'lib_strategy',
                    'lib_source'
                ]) + "\n"
            )
        sys.stdout.flush()

        for tmpdict in results:
            sys.stdout.write(
                "\t".join([
                        tmpdict['runacc']    ,
                        tmpdict['sciname']   ,
                        tmpdict['platname']  ,
                        tmpdict['ispublic']  ,
                        tmpdict['strategy']  ,
                        tmpdict['lib_source']
                    ]) + "\n"
            )
            sys.stdout.flush()

def main():
    opts    = getOpts()

    countersamples = []

    if opts.counterspps:        
        with open(opts.counterspps, 'r') as f:
            for i in f.readlines():
                countersamples.append(i.strip())
                
    myclass = entrez(term = opts.term,
                     db   = "sra",
                     type = "docsum")

    mytype  = myclass._get_type()

    if not mytype:
        sys.stdout.write("\nNo data for %s\n" % opts.term)
        sys.stdout.flush()
        exit()

    results = formatmytype(mytype=mytype, countersamples=countersamples)

    printthisout(results)

if __name__ == "__main__":
    main()



