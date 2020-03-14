# geneGet

## Installation

Using `git`:

```Shell
git clone https://github.com/Ulises-Rosas/geneTable.git
cd geneTable
python3 setup.py install
```

## Usage

```Bash
getFeatures "Alopias vulpinus" -p
```
![](https://github.com/Ulises-Rosas/geneTable/blob/master/img/Alopias_vulpinus_getFeatures.png)

Filnames are composed by using species name as well as its title by default. Horizontal line depicts three sequences. If there were more than 200 hundred NCBI ids per species, the number of downloaded tables per species is controled with the argument `--cache`.
