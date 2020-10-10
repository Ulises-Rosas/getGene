
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
                 request = None,
                 threads = 1):

        self.request = request
        self.threads = threads
        # self.out_dir  = out_dir
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

        # accession  = 'GCF_001319825.1'
        # spps = 'Yersinia wautersii'
        # mydirname = (spps + "__" + accession).replace(" ", "_")

        down_url = "/".join([
                         self.base_url, 
                         self.ass_acc, 
                         self.data  % accession
                        ])

        sys.stdout.write("Downloading: %s, %s\n" % (spps, accession))
        sys.stdout.flush()

        response = requests.get( down_url, headers = self.header ) 

        start = time.time()
        while not response.ok:
            time.sleep(0.5)

            final = time.time()
            if (final - start) > timeout:
                sys.stderr.write("Timeout response: %s, %s\n" % (spps, accession))
                sys.stderr.flush()
                return None
        
        myzip = zipfile.ZipFile( io.BytesIO(response.content) )
        myzip.extractall(self.sppsdir)

        return 1

    def iterate_genome(self):

        # Multiprocess for moving files
        with Pool(processes = self.threads) as p:

            for acc, spps in self.request:
                # acc = 'GCF_001319825.1'
                # spps = 'Yersinia wautersii'
                self.sppsdir = (spps + "__" + acc).replace(" ", "_")
                stat_ok      = self._get_genome(accession=acc, spps=spps)

                if not stat_ok:
                    continue                

                myfiles  = glob.glob(
                                os.path.join(
                                    self.sppsdir,
                                    "ncbi_dataset",
                                    "data",
                                    acc,
                                    "*_genomic.fna"))

                [ *p.map(self._move_file, myfiles) ]

                shutil.rmtree(
                        os.path.join(
                            self.sppsdir,
                            "ncbi_dataset"))

                readmefil = os.path.join(self.sppsdir, 'README.md')

                if os.path.exists(readmefil):
                    os.remove(readmefil)
                
    def taxon_descriptor(self,
                         timeout = 60,
                         spps    = None):
        # spps    = 'Yersinia wautersii'
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
            title    = checkkey(org, 'title')

            out.append({
                'accession'  : accession,
                'seq_length' : seq_length,
                'conting_n50': n50,
                'ass_level'  : ass_level,
                'sci_name'   : sci_name,
                'title'      : title
                })

        return out

