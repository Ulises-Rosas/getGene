#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 12:48:35 2020

@author: ulisesrosasp
"""

import re
import sys
from collections import Counter

getMaxNu = lambda d: sorted([len(v) for _,v in d.items()], reverse = True)[0]

def genomeDS(docsum):
#    docsum = rowString 
    colpat   = '^.+Name="(.+)" .+>.{0,}</Item>$'
    valpat   = '^.+Name=".+" .+>(.{0,})</Item>$'

    fulfill  = lambda d,n: {k: v + ( [""] * (n - len(v)) ) for k,v in d.items() }

    sdocsum  = docsum.split("\n")
    
    out = {}
    for i in sdocsum:

        if not re.findall("</Item>$", i):
            continue

        colname = re.sub(colpat, "\\1", i)
        val     = re.sub(valpat, "\\1", i)

        out[colname] = out[colname] + [val] if out.__contains__(colname) else [val]
        
    return fulfill( out, getMaxNu(out) )

def dictToPrint(dic = None, sep = "\t"):
    """
    associated to getGenomes script
    """
#    dic = genomeDS(rawDoc)
    head = list(dic)
    
    out  = [sep.join(head)]
    nrow = len(dic[head[0]])
    
    for p in range(0, nrow):
        row = []
        for h in head:
            row += [dic[h][p]]
            
        out += [sep.join(row)]
        
    return out

def dictToCSV(dic = None, outname = None, sep = ",", quiet = False):
    # c = {'coi': 149, 'hsp70': 7, 'trna-ser': 6, 'trna-leu': 6, '12s ribosomal rna': 4, 'cytb': 3, '16s ribosomal rna': 3, 'trna-pro': 3, 'trna-ile': 3, 'trna-glu': 3, 'trna-his': 3, 'trna-asp': 3, 'trna-trp': 3, 'trna-val': 3, 'trna-gly': 3, 'trna-phe': 3, 'trna-thr': 3, 'trna-cys': 3, 'trna-asn': 3, 'trna-ala': 3}
    # dic = c
    outputname = outname + ".txt"

    f = open(outputname, "w")
    for k,v in dic.items():
        f.write("{gene}{sep}{freq}\n".format(gene = k, sep = sep, freq = v))
    f.close()

    if not quiet:
        sys.stdout.write("\n")
        sys.stdout.write("Table stored at %s\n" % outputname)


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
            
