#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 12:48:35 2020

@author: ulisesrosasp
"""

import re


getMaxNu = lambda d: sorted([len(v) for _,v in d.items()], reverse = True)[0]

def genomeDS(docsum):
#    docsum = rawDoc
    
    pat1     = '^.+Name="(.+)" .+>(.{0,})</Item>$'
    
    nameItem = lambda s: re.sub(pat1, "\\1#\\2", s) if re.findall("</Item>$", s) else None
    fulfill  = lambda d,n: {k: v + ( [""] * (n - len(v)) ) for k,v in d.items() }
        
    sdocsum  = docsum.split("\n")
    
    out = {}
    for i in filter(None, map(nameItem, sdocsum)):
        
        colname,val  = i.split("#")                        
        out[colname] = out[colname] + [val] if out.__contains__(colname) else [val]
        
    return fulfill( out, getMaxNu(out) )

def dictToPrint(dic):
#    dic = genomeDS(rawDoc)
    head = list(dic)
    
    out  = ["\t".join(head)]
    nrow = len(dic[head[0]])
    
    for p in range(0, nrow):
        row = []
        for h in head:
            row += [dic[h][p]]
            
        out += ["\t".join(row)]
        
    return out
            
