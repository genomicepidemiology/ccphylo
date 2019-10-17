# Getting Started #

```
git clone https://bitbucket.org/genomicepidemiology/ccphylo.git
cd ccphylo && make

./ccphylo -h
./ccphylo dist -i path/to/kma/alignments/*.mat.gz -o results.phy
./ccphylo tree -i results.phy -o results.nwck
./ccphylo union -i path/to/kma/alignments/*.res
```

# Introduction #
CCPhylo is a software suite designed to easen phylogenitic analysis based on alignments from KMA. 
Neighbour-Joining can be run without KMA alignments, if a Phylip distance matrix is given using "ccphylo tree".

If you use CCPhylo for your published research, then please cite:
*Article on the way*


# Usage #
For practical reasons you might want to add ccphylo to your path, this is usually done with:

```
mv ccphylo ~/bin/
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
1. *Article on the way*

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
