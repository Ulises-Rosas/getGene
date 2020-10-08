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

### getgenomes

Look and download genomes

```Shell
getgenomes Stomiatiformes
```
```
Organism_Name   Assembly_Accession
Borostomias antarcticus GCA_900323325.1
```

### looksra

Look SRA information

```Shell
looksra Argentiniformes
```
```
accession       species_name    platform        is_public       lib_strategy    lib_source
SRR11679483     Argentina silus Illumina HiSeq 4000     true    WGS     GENOMIC
SRR11537143     Argentina sphyraena     Illumina HiSeq 4000     true    WGS     GENOMIC
ERR3332509      Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332508      Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332507      Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
ERR3332506      Opisthoproctus soleatus Illumina HiSeq X Ten    true    WGS     GENOMIC
SRR5997680      Argentina sp. CUR14063.G        Illumina HiSeq 2000     true    RNA-Seq TRANSCRIPTOMIC
```
