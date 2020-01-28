

import re
import urllib.request
import argparse
import matplotlib.pyplot as plt
from collections import Counter

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

        self.efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db="

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

    def get_seqs(self,
                 ids=""):

        if ids == "":
            ids = self._get_ids()

        # delete this
        # self = entrez(term= "JL554673:JL597291", termType= "[ACCN]", db = "nuccore", type = "fasta")
        # ids0 = self._get_ids()
        # ids = ids0[0:4]

        i = 0
        if self.printing:
            while (i <= len(ids)):
                complete_efetch_url = self.efetch_url + self.db +\
                                      "&id=" + ",".join( ids[i:i + self.cache] ) +\
                                      "&rettype=" + self.type

                print(urllib.request.urlopen(complete_efetch_url).read().decode('utf-8'))
                i += self.cache
        else:
            string = ""
            while (i <= len(ids)):
                complete_efetch_url = self.efetch_url + self.db +\
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

        def GenesftByText(page, keyWords=["gene"]):

            feats = list(filter(None, page.split(">Feature ")))
            # keyWords = ["gene", "rRNA", "tRNA"]

            matchPattern = "[0-9<>\t]+%s\n[\t]+[A-Za-z]+\t.+?(?=\n)"
            subPattern = "([0-9<>]+)\t([0-9<>]+)\t%s\n[\t]+[A-Za-z]+\t(.*)"

            allFeats = []

            for ft in feats:
                # ft = feats[8]
                # keyWords = ["gene", "rRNA"]
                tmpRegions = []

                for key in keyWords:

                    tmpMatch = re.findall(matchPattern % key, ft)

                    for mtchs in tmpMatch:
                        tmpSub = re.sub(subPattern % key
                                        , "\\1,\\2,\\3,%s" % key
                                        , mtchs).replace("<", "").replace(">", "")

                        tmpRegions.append(tmpSub)

                positions = [",".join(i.split(',')[0:2]) for i in tmpRegions]

                for josp in list(set(positions)):
                    josr = [i for i in tmpRegions if re.findall(josp, i)]

                    lenOfRegionName = [len(i.split(',')[2]) for i in josr]

                    shortestWord = [x for _, x in sorted(zip(lenOfRegionName, josr))][0]

                    allFeats.append(

                        shortestWord.split(',')[2].lower()
                    )

            return dict(Counter(allFeats))

        #dict1 = {}
        superPage = ''

        if ids0.__len__() <= 200:

            complete_efetch_url = self.efetch_url + \
                                  self.db + \
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
                                      self.db + \
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

            return dict(sortedDict1) if cutOff == None else dict(sortedDict1[0:int(cutOff)])

