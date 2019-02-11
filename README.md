# geneTable
Gene availability of given species by using feature table from NCBI's E-utilities

## Requirements
* Python >=3.5
* Matplotlib

## Bar chart of gene availabilty 

```Bash
python3 ./geneTable.py -type "ft" --plot "Oreochromis mossambicus" --cutOff 20
```
![](https://github.com/Ulises-Rosas/geneTable/blob/master/img/Oreochromis_mossambicus_GeneAvailability.png)

Filnames are composed by using species name as well as its title by default. If there were more than 200 hundred NCBI ids per species, the number of downloaded tables per species is controled with the argument `-cache`
