[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ccphylo/README.html)
# Getting Started #

```
git clone https://bitbucket.org/genomicepidemiology/ccphylo.git
cd ccphylo && make

./ccphylo -h
./ccphylo tree test.phy.gz
```

# Introduction #
CCPhylo is a software suite designed to easen phylogenetic analysis based on alignments from KMA. 
Neighbor-Joining can be run without KMA alignments, if a Phylip distance matrix is given using "ccphylo tree".

If you use CCPhylo for your published research, then please cite:
Philip T.L.C. Clausen
"Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining",
Bioinformatics, 2023.
Malte B. Hallgren, Søren Overballe-Petersen, Ole Lund, Henrik Hasman & Philip T.L.C. Clausen
"MINTyper: an outbreak-detection method for accurate and rapid SNP typing of clonal clusters with noisy long reads",
Biology, Methods & Protocols 2021.


# Usage #

```
./ccphylo -h
./ccphylo dist -i path/to/kma/alignments/*.mat.gz -o results.phy
./ccphylo tree -i results.phy -o results.nwck
./ccphylo union -i path/to/kma/alignments/*.res
```

For practical reasons you might want to add ccphylo to your path, this is usually done with:

```
mv ccphylo ~/bin/
```

# Examples #
CCPhylo tree has a variety of options to construct hierarchical clusterings in newick format, including the Dynamic 
Neighbor-Joining (dnj, default) algorithm that produces exact Neighbor-Joining trees in near quadric time, other 
methods may be applied using the "-m" option. Both relaxed and strict Phylip format are accepted in both lower 
triangular and full matrix form, all is auto detected (please note that tabs are not allowed in the entry names, 
while spaces are accepted).
Example of creating dnj tree, where "distances.phy" is the phylip distance matrix.
```
ccphylo tree -i distances.phy
```

Generate a distance matrix and a tree between the samples, using "CP015990.1" as a reference, and include a matrix 
describing how many nucleotide positions that were used between each sample.
```
ccphylo dist -i path/to/kma/alignments/*.mat.gz -o path/to/ccphylo/results/result.phy -n path/to/ccphylo/results/result.num -r CP015990.1
ccphylo tree -i path/to/ccphylo/results/result.phy -o path/to/ccphylo/results/result.nwck
```

Find all pairwise shared templates between samples, and make a diatance matrix and tree for each of them.
```
ccphylo union -i path/to/kma/alignments/*.res -o path/to/ccphylo/results/result.union -B path/to/kma/db
ccphylo dist -i path/to/ccphylo/results/result.union -o path/to/ccphylo/results/result.phy
ccphylo tree -i path/to/ccphylo/results/result.phy -o path/to/ccphylo/results/result.nwck
```

This can also be done without saving the intermediate results of CCPhylo, with:
```
ccphylo union -i path/to/kma/alignments/*.res -B path/to/kma/db | ccphylo dist | ccphylo tree
```


# Installation Requirements #
In order to install CCPhylo, you need to have a C-compiler and zlib development files installed.
Zlib development files can be installed on unix systems with:
```
sudo apt-get install libz-dev
```

# Help #
Usage and options are available with the "-h" option on all programs in the software suite.
If in doubt, please mail any concerns or problems to: *plan@dtu.dk*.

# Citation #
1. Philip T.L.C. Clausen "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining", Bioinformatics, 2023.
2. Malte B. Hallgren, Søren Overballe-Petersen, Ole Lund, Henrik Hasman & Philip T.L.C. Clausen "MINTyper: an outbreak-detection method for accurate and rapid SNP typing of clonal clusters with noisy long reads", Biology, Methods & Protocols 2021.

# License #
Copyright (c) 2017, Philip Clausen, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
