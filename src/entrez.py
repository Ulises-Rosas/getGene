

import re
import sys
import time
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
        self.preterm = term.replace(" ", "%20")
        # self.term = term.replace(" ", "%20") 
        self.gene_string = gene_string

        self.db = db
        self.cache = cache
        self.printing = printing
        self.Lmin = Lmin
        self.Lmax = Lmax

        self.efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=" + self.db

        self.ids = []

    @property
    def term(self):

        if self.db != 'nuccore':
            return self.preterm

        termType = "[Organism]" if not re.findall("\\[", self.preterm) else ""

        if self.gene_string:
            genes       = self.gene_string.split(",")
            genebool    = " OR ".join([i + "[All Fields]" for i in genes])
            gene_string = genebool if len(genes) == 1 else "("+genebool+")"

        else:
            gene_string = self.gene_string

        if self.Lmin != "" and self.Lmax != "":
            Lrange = "(" + str(self.Lmin) + "[SLEN] :" + str(self.Lmax) + "[SLEN])"
        else:
            Lrange = ""

        myopts = filter(None, [self.preterm + termType, gene_string, Lrange])

        return " AND ".join(myopts).replace(" ", "%20")

    @property
    def esarch_url(self):

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="
        return url + self.db + "&term=" + self.term

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

        self.ids = [re.sub("<Id>([0-9\\.]+)</Id>", "\\1", i) for i in re.findall("<Id>[0-9\\.]+</Id>", ids_page)]

        return self.ids

    def _get_type(self, type = None):

        if type is not None:
            self.type = type

        ids = self._get_ids() if not self.ids else self.ids

        if not ids:
            return None
            
        out = []
        i   = 0
        totallength = len(ids)

        while i <= totallength:

            if totallength > self.cache:
                sys.stderr.write("\rDownloading metadata (%6.2f%%)" % (i*100/totallength) )
                sys.stderr.flush()

            myids      = ",".join( ids[i:i + self.cache] )
            target_url = self.efetch_url + "&id=" + myids + "&rettype=" + self.type
            tmp_out    = urllib.request.urlopen(target_url).read().decode('utf-8')
            time.sleep(0.75)

            out += [tmp_out]
            i   += self.cache

        if totallength > self.cache:
            sys.stderr.write('\n')
            sys.stderr.flush()
            
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

    def getMaxNu(self, mydict):

        return sorted([len(v) for _,v in mydict.items()], reverse = True)[0]

    def genomeDS(self):

        docsum = list(filter(None,self._get_type()))

        if not docsum:
            return None

        colpat   = '^.+Name="(.+)" .+>.{0,}</Item>$'
        valpat   = '^.+Name=".+" .+>(.{0,})</Item>$'

        fulfill  = lambda d,n: {k: v + ( [""] * (n - len(v)) ) for k,v in d.items() }

        # superout = {}
        out = {}

        for ds in docsum:

            sdocsum  = ds.split("\n")
                
            for i in sdocsum:

                if not re.findall("</Item>$", i):
                    continue

                colname = re.sub(colpat, "\\1", i)
                val     = re.sub(valpat, "\\1", i)

                out[colname] = out[colname] + [val] if out.__contains__(colname) else [val]
            

        return fulfill(out, self.getMaxNu(out)) 

        # return superout


# self = entrez(term = 'Clupeiformes',
#                     db   = "genome",
#                     type = "docsum",
#                     cache=200)
# self.type = 'docsum'