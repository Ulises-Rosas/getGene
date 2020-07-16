
import re
import os
import io
import sys
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

        down_url = os.path.join( self.base_url, 
                                 self.ass_acc, 
                                 self.data  % accession )

        response = requests.get( down_url, headers = self.header ) 

        start = time.time()
        while not response.ok:
            time.sleep(0.5)

            final = time.time()
            if (final - start) > timeout:
                sys.stdout.write("\nTimeout response: %s" % spps)
                sys.stdout.flush()
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
