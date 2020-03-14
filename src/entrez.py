

import re
import argparse
import urllib.request
from collections import Counter
import matplotlib.pyplot as plt
from genetable.utils import GenesftByText

class entrez:

    def __init__(self,
                 term="",
                 type="",
                 db="",
                 gene_string="",
                 Lmin="",
                 Lmax="",
                 cache = 200,
                 printing = True):

        self.type = type
        self.term = term.replace(" ", "%20") 
        self.db = db
        self.cache = cache
        self.printing = printing

        termType = "[Organism]" if not re.findall("\[", self.term) else ""

        if self.db == "nuccore":

            if gene_string != "" and len(gene_string.split(",")) == 1:
                gene_string = " OR ".join([i + "[All Fields]" for i in gene_string.split(",")])

            elif gene_string != "" and len(gene_string.split(",")) > 1:
                gene_string = "(" + \
                       " OR ".join([i + "[All Fields]" for i in gene_string.split(",")]) + \
                       ")"

            if Lmin != "" and Lmax != "":
                Lrange = "(" + str(Lmin) + "[SLEN] :" + str(Lmax) + "[SLEN])"
            else:
                Lrange = ""

            self.term = re.sub(" ",
                               "%20",
                               " AND ".join(
                                   [i for i in [self.term + termType, gene_string, Lrange] if i != ""]
                                    )
                               )

        self.esarch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=" + \
                          self.db + "&term=" + self.term

        self.efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=" + self.db

        self.ids = []

    def _get_ids(self):

        page = None

        while page is None:
            try:
                page = urllib.request.urlopen(self.esarch_url).read().decode('utf-8')

            except urllib.error.HTTPError:
                pass

        counts = re.sub(".*><Count>([0-9]+)</Count><.*", "\\1", page.replace("\n", ""))

        complete_esearch_url = self.esarch_url + "&retmax=" + counts

        ids_page = urllib.request.urlopen(complete_esearch_url).read().decode('utf-8')

        self.ids = [re.sub("<Id>([0-9\.]+)</Id>", "\\1", i) for i in re.findall("<Id>[0-9\.]+</Id>", ids_page)]

        return self.ids

    def _get_type(self, type = None):

        if type is not None:
            self.type = type

        ids = self._get_ids() if not self.ids else self.ids

        if not ids:
            return None
            
        out = ""
        i   = 0

        while i <= len(ids):

            target_url = self.efetch_url + \
                         "&id=" + ",".join( ids[i:i + self.cache] ) + \
                         "&rettype=" + self.type
    
            tmp_out    = urllib.request.urlopen(target_url).read().decode('utf-8')

            out += tmp_out
            i   += self.cache

        return out

    def get_seqs(self,
                 ids=""):

        if ids == "":
            ids = self._get_ids()

        # delete this
        # ids0 = self._get_ids()
        # ids = ids0[0:4]

        i = 0
        if self.printing:
            while (i <= len(ids)):
                complete_efetch_url = self.efetch_url +\
                                      "&id=" + ",".join( ids[i:i + self.cache] ) +\
                                      "&rettype=" + self.type

                print(urllib.request.urlopen(complete_efetch_url).read().decode('utf-8'))
                i += self.cache
        else:
            string = ""
            while (i <= len(ids)):
                complete_efetch_url = self.efetch_url +\
                                      "&id=" + ",".join(ids[i:i + self.cache]) +\
                                      "&rettype=" + self.type

                page = urllib.request.urlopen(complete_efetch_url).read().decode('utf-8')
                string += page
                i += self.cache

            return string

    def feature_table(self, keyword, cutOff):

        #keyword = ["gene"]
        #self = entrez(term= "Litopenaeus vannamei", db= "nuccore", type= "ft")
        #self = entrez(term="Anguilla anguilla", db="nuccore", type="ft")

        ids0 = self._get_ids()
        #ids.__len__()
        #ids0 = ids[0:2000]

        #dict1 = {}
        superPage = ''

        if ids0.__len__() <= 200:

            complete_efetch_url = self.efetch_url + \
                                  "&id=" + \
                                  ",".join(ids0) + \
                                  "&rettype=" + \
                                  self.type
            page = None

            while page is None:
                try:
                    page = urllib.request.urlopen(complete_efetch_url).read().decode('utf-8')

                except urllib.error.HTTPError:
                    pass

            superPage += page

        else:

            i = 0
            while(  ids0.__len__() > i ):

                complete_efetch_url = self.efetch_url + \
                                      "&id=" + \
                                      ",".join(ids0[i:i + self.cache]) + \
                                      "&rettype=" + \
                                      self.type
                page = None

                while page is None:
                    try:
                        page = urllib.request.urlopen(complete_efetch_url).read().decode('utf-8')

                    except urllib.error.HTTPError:
                        pass

                superPage += page

                i += self.cache

                if i > ids0.__len__():

                    print("progress: {0}/{0} tables downloaded".format( ids0.__len__() ) )
                else:

                    print("progress: {0}/{1} tables downloaded".format(i, ids0.__len__()))


        dict1 = GenesftByText(page=superPage,
                              keyWords=keyword)

        if dict1.__len__() == 0:

            return None
        else:
            sortedDict1 = sorted(dict1.items()
                                 , key     = lambda kv: kv[1]
                                 , reverse = True)

            return dict(sortedDict1) if cutOff is None else dict(sortedDict1[0:int(cutOff)])

