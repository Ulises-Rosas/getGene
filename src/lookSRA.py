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

def checkmatch(mystr, exists, extract):

    return re.sub(extract, "\\1", mystr) if re.findall(exists,mystr) else ""

def formatmytype(mytypes, countersamples):
    out = []

    for mytype in mytypes:

        tree = ET.fromstring(mytype)

        for child in tree:

            expxml = child[1].text
            runa   = child[2].text

            if expxml:
                # libsrcpatt = '.*<LIBRARY_SOURCE>(.+)</LIBRARY_SOURCE>.*'
                # srcname    = re.sub(libsrcpatt, "\\1", expxml)
                # print(expxml)
                srcname = checkmatch(
                            expxml, 
                            '<LIBRARY_SOURCE>',
                            '.*<LIBRARY_SOURCE>(.+)</LIBRARY_SOURCE>.*'
                        )
                # libstrpatt = '.*<LIBRARY_STRATEGY>(.+)</LIBRARY_STRATEGY>.*'
                # straname = re.sub(libstrpatt, "\\1", expxml)
                straname = checkmatch(
                                expxml,
                                '<LIBRARY_STRATEGY>',
                                '.*<LIBRARY_STRATEGY>(.+)</LIBRARY_STRATEGY>.*'
                            )

                rawlayout = checkmatch(
                                expxml,
                                '<LIBRARY_LAYOUT>',
                                '.*<LIBRARY_LAYOUT>(.+)</LIBRARY_LAYOUT>.*'
                            )
                layout = re.sub("[ ]{0,1}<([A-Za-z]+).*","\\1",rawlayout).lower()
                
                
                # platpatt = '.*instrument_model="(.*?)".*'
                # platname = re.sub(platpatt, "\\1", expxml)
                platname = checkmatch(
                                expxml,
                                'instrument_model=',
                                '.*instrument_model="(.*?)".*'
                            )
                # scipatt = '.*ScientificName="(.*?)".*'
                # sciname = re.sub(scipatt, "\\1", expxml)
                sciname = checkmatch(
                                expxml,
                                'ScientificName=',
                                '.*ScientificName="(.*?)".*'
                            )

                if countersamples:
                    if (sciname in countersamples) or (not sciname):
                        expxml = None

            if runa:
                # accpatt = '.* acc="(.*?)".*'
                # runacc  = re.sub(accpatt, "\\1", runa)
                runacc  =  checkmatch(
                                runa,
                                'acc=',
                                '.* acc="(.*?)".*'
                            )

                # ispupatt = '.* is_public="(.*?)".*'
                # ispublic = re.sub(ispupatt, "\\1", runa)
                ispublic = checkmatch(
                                runa,
                                'is_public=',
                                '.* is_public="(.*?)".*'
                            )

                # basepatt = '.* total_bases="(.*?)".*'
                # totalbase = re.sub(basepatt, "\\1", runa)
                totalbase = checkmatch(
                                runa,
                                'total_bases=',
                                '.* total_bases="(.*?)".*'
                            )
            
            if expxml and runa:
                tmp = {
                    'sciname'    : sciname,
                    'platname'   : platname,
                    'runacc'     : runacc,
                    'ispublic'   : ispublic,
                    'strategy'   : straname,
                    'lib_source' : srcname,
                    'total_bases': totalbase,
                    'layout'     : layout
                }
                out.append(tmp)

    return out

def printthisout(results):

    if results:
        sys.stdout.write(
            "\t".join([
                    'accession'   ,
                    'total_bases' ,
                    'layout'      ,
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
                        tmpdict['runacc']     ,
                        tmpdict['total_bases'],
                        tmpdict['layout']     ,
                        tmpdict['sciname']    ,
                        tmpdict['platname']   ,
                        tmpdict['ispublic']   ,
                        tmpdict['strategy']   ,
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

    results = formatmytype(mytypes=mytype, countersamples=countersamples)

    printthisout(results)

if __name__ == "__main__":
    main()



