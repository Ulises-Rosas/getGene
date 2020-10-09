
import re
import os
import io
import sys
import json
import time
import zipfile
import requests
from urllib import parse

# from genetable.entrez import entrez
# from genetable.utils import genomeDS, dictToPrint


import glob
import shutil
from multiprocessing import Pool

class Datasets:
    def __init__(self,
                 dictstring = None,
                 out_dir = None,
                 threads = 1):

        self.dictstring = dictstring
        self.threads    = threads
        self.out_dir    = out_dir
        # self.taxon      = taxon

        self.base_url = "http://api.ncbi.nlm.nih.gov/datasets/v1alpha"
        self.ass_acc  = "download/assembly_accession"
        self.data     = "%s?include_sequence=true&hydrated=FULLY_HYDRATED"
        self.header   = {"accept": "application/zip"}


        # placeholder
        self.sppsdir = ""

    def _move_file(self, file):
        shutil.copy(file, self.sppsdir)

    def _get_genome(self,
                    timeout   = 60,
                    accession = None,
                    spps      = None):

        down_url = "/".join([
                         self.base_url, 
                         self.ass_acc, 
                         self.data  % accession
                        ])

        response = requests.get( down_url, headers = self.header ) 

        start = time.time()
        while not response.ok:
            time.sleep(0.5)

            final = time.time()
            if (final - start) > timeout:
                sys.stderr.write("\nTimeout response: %s" % spps)
                sys.stderr.flush()
                return

        sys.stdout.write("\nDownloaded genome: %s" % spps)
        sys.stdout.flush()
        
        myzip = zipfile.ZipFile( io.BytesIO(response.content) )
        myzip.extractall()

    def iterate_genome(self):

        mydict = self.dictstring
        out_dir = self.out_dir

        if not os.path.isdir( out_dir ):
            os.mkdir( out_dir )

        with Pool(processes = self.threads) as p:

            for n, spps in enumerate( mydict['Organism_Name'] ):

                ass_acc = mydict['Assembly_Accession'][n]

                if not ass_acc:
                    continue

                self._get_genome(accession=ass_acc, spps=spps)
                self.sppsdir = os.path.join( out_dir, spps.replace(" ", "_") )

                if not os.path.isdir( self.sppsdir ):
                    os.mkdir( self.sppsdir )

                myfiles  = glob.glob(
                            os.path.join(
                                "ncbi_dataset",
                                "data",
                                ass_acc,
                                "*_genomic.fna"
                            )
                        )

                [ *p.map(self._move_file, myfiles) ]

                shutil.rmtree("ncbi_dataset")

    def taxon_descriptor(self,
                         timeout = 60,
                         spps    = None):
        # spps    = 'Argentiniformes'
        req_url = "/".join([
                        self.base_url,
                        'genome/taxon/%s?limit=all&returned_content=COMPLETE' % spps.replace(" ", "%20")
                        ])

        response = requests.get( req_url, headers = {"accept": "application/json"} ) 

        start = time.time()
        while not response.ok:
            time.sleep(0.5)

            tmpload = {}
            try:
                tmpload = json.loads(response.content)
            except json.JSONDecodeError:
                continue

            if tmpload:
                if tmpload.__contains__('error'):
                    sys.stderr.write(tmpload['message'] + "\n")
                    sys.stderr.flush()
                    return None
                else:
                    continue

            final = time.time()
            if (final - start) > timeout:
                sys.stderr.write("\nTimeout response: %s" % spps)
                sys.stderr.flush()
                return None
        # sys.stdout.write("\nDownloaded genome: %s" % spps)
        # sys.stdout.flush()
        loaded = json.loads(response.content)

        if not loaded:
            sys.stderr.write("\nNo data for %s\n" % spps)
            sys.stderr.flush()
            return None

        out = []

        checkkey = lambda json, key: json[key] if json.__contains__(key) else ''

        for li in loaded['assemblies']:

            assem = checkkey(li,'assembly')

            if not assem:
                continue

            accession  = checkkey(assem, 'assembly_accession')
            seq_length = checkkey(assem, 'seq_length')
            n50        = checkkey(assem, 'contig_n50')
            ass_level  = checkkey(assem, 'assembly_level')

            org = checkkey(assem, 'org')

            if not org:
                continue
            sci_name = checkkey(org, 'sci_name')

            out.append({
                'accession'  : accession,
                'seq_length' : seq_length,
                'conting_n50': n50,
                'ass_level'  : ass_level,
                'sci_name'   : sci_name
                })

        return out


# self = Datasets(out_dir="genome_out")




# if __name__ == "__main__":
    # term = "yersinia"
    # rawOut = entrez(
    #             term = term, # opts.term
    #             db   = "genome" ,
    #             type = "docsum")

    # rawString = rawOut._get_type()

    # if rawString is not None:
    #     # rowString  = rawOut._get_type()
    #     dictString = genomeDS(rawString)

    #     Datasets(dictString , out_dir= "genome_out").iterate_genome()
    # else:
    #     sys.stdout.write("\nCheck term: %s" % term)
    #     sys.stdout.flush()

    # sys.stdout.write("\n")
