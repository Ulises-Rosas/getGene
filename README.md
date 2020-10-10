# geneGet

## Installation

Using `git`:

```Shell
git clone https://github.com/Ulises-Rosas/geneTable.git
cd geneTable
python3 setup.py install
```

## Usage


### getFeatures

Get gene availability and plot it (i.e. `-p`).

```Shell
getFeatures "Alopias vulpinus" -p
```
![](https://github.com/Ulises-Rosas/geneTable/blob/master/img/Alopias_vulpinus_getFeatures.png)

Filnames are composed by using species name as well as its title by default. Horizontal line depicts three sequences. If there were more than 200 hundred NCBI ids per species, the number of downloaded tables per species is controled with the argument `--cache`.

### lookgenomes

Look at genome metadata

```Shell
lookgenomes Clupeiformes
```
```
accession       seq_length      conting_n50     sci_name        ass_level
GCA_007927625.1 812159898       1632759 Coilia nasus    Chromosome
GCF_900700415.1 725670187       1151065 Clupea harengus Chromosome
GCA_900700415.1 725670187       1151065 Clupea harengus Chromosome
GCA_902175115.3 790426535       22496   Clupea harengus Scaffold
GCA_900323705.1 789621745       28410   Clupea harengus Scaffold
GCA_000966335.1 807695261       22277   Clupea harengus Scaffold
GCA_900499035.1 949617276       9398    Sardina pilchardus      Scaffold
GCA_003604335.1 641491539       10878   Sardina pilchardus      Scaffold
GCA_003651195.1 815647530       129889  Tenualosa ilisha        Scaffold
GCA_004329175.1 710279582       1092    Tenualosa ilisha        Scaffold
GCF_900700375.1 567401054       3059612 Denticeps clupeoides    Chromosome
GCA_900700375.2 567401054       3059612 Denticeps clupeoides    Chromosome
GCA_900700375.1 567401054       3059612 Denticeps clupeoides    Chromosome
GCA_900700345.2 457704509       384267  Denticeps clupeoides    Scaffold
```

### looksra

Look at SRA metadata

```Shell
looksra Argentiniformes
```
```
accession       total_bases     layout  species_name    platform        is_public       lib_strategy    lib_source
SRR11679483     1676957378      paired  Argentina silus Illumina HiSeq 4000     true    WGS     GENOMIC
SRR11537143     2577408128      paired  Argentina sphyraena     Illumina HiSeq 4000     true    WGS     GENOMIC
ERR3332509      26274979990     paired  Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332508      26507370802     paired  Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332507      32283025068     paired  Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332506      25164048528     paired  Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
SRR5997680      4591698120      paired  Argentina sp. CUR14063.G        Illumina HiSeq 2000     true    RNA-Seq TRANSCRIPTOMIC
```

### getgenomes

Get genomes with a file containing accession codes and species names. For example, given the following file `listacc.txt` :

```
cat listacc.txt
```
```
GCF_001319825.1,Yersinia wautersii
GCF_000493475.1,Yersinia wautersii
```
We can get genomes from that file by using:

```
getgenomes listAcc.txt 

head Yersinia_wautersii__GCF_001319825.1/GCF_001319825.1_5139_1_1_genomic.fna
```
```
>NZ_CVMG01000001.1 Yersinia wautersii strain WP-931201, whole genome shotgun sequence
ATGGCCCACGGTGGAAAACTGGCCACAGGTTGAGTTGCCGGAACTGCCGCAATGGCTTTTGATTGAAGCGGTCAATCAGG
GTTATATTGTTCCCGACTGGCCGCCAGTCGTATAGGCCTGCCCAACAACCCCTCCAGTCGGGGTTGTTGGTTTCTCTGTT
GTGCCACCCCCCACACAATCCTCATCACCTGCCCCGCGCGCAGTAATCCGGCATCATAGCGAATGAACGCTTAACCGGAG
AAAAACGCATGTCTGCAACCGATTACCACCACGGTGTGCGCGTCATTGAAATCAGCGAAGGCACTCGCCCGATCCGCACT
GTCAGTACGGCGGTAGTCGGGATGGTCTGTACTTCCGATGATGCTGACCCCACTCTGTTCCCACTCAATACCCCGGTATT
ACTCACCGATGTGCTGGCCGCCAGCGGCAAGGCCGGTGAAACCGGCACATTAGCCCATTCACTGGATGCTATCAGCGACC
AAACCAAGCCCGTGACTATTGTTGTCCGCGTGGCTCAGGGGGAAACCGAAGCCGAAACTACCTCCAATATTATCGGCGGC
TCCACGCCAGATGGCCGTTATACCGGCATGAAAGCGCTGTTAGCGGCGCAGGGTAAGTTTGACGTCAAGCCCCGTATTTT
AGGGGTGCCCGGTCATGACACTCTGGCGGTATCCACTGAGCTACTTTCCATCGCTCAGAGCCTACGTGCCTTTGCCTACA
```

