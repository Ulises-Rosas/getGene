# geneTable
Gene availability of given species by using feature table from NCBI's E-utilities

## Requirements
* Python >=3.5
* Matplotlib

#### Installation

Using `git`:

```Shell
git clone https://github.com/Ulises-Rosas/geneTable.git
cd geneTable
python3 setup.py install
```

## Bar chart of gene availabilty 

```Bash
geneTable.py -type "ft" --plot "Oreochromis mossambicus" --cutOff 20
```
![](https://github.com/Ulises-Rosas/geneTable/blob/master/img/Oreochromis_mossambicus_GeneAvailability.png)

Filnames are composed by using species name as well as its title by default. Horizontal line depicts three sequences. If there were more than 200 hundred NCBI ids per species, the number of downloaded tables per species is controled with the argument `-cache`
