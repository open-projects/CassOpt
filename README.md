### The project was supported by RFBR (the RFBR project number is 17-04-02186).

# CassOpt (Cassette Optimizer)
The program creates an optimized mini-gene translated sequence using a set of sub-sequences of predicted Minor Histocompatibility Antigenes (MiHA). The optimized mini-genes are used for experemental bulk validation of immunogenicity of predicted MiHAs. The translated MiHA coding gene regions with flanks are arranged to a cassette in the order where there are no immunogenic flank junction sequences. It shortens the translated mini-gene sequence and prevents false positive results of the experiments.

## Requirements
* Linux or MacOS
* python >= 3.5
* netMHCpan >= 4.0
* tcsh >= 6.18 (required for netMHCpan)

## Installation
* install netMHCpan4 program which was described in the paper:

  J Immunol. 2017 Nov 1;199(9):3360-3368. doi: 10.4049/jimmunol.1700893. Epub 2017 Oct 4.
  
  (https://www.jimmunol.org/content/early/2017/10/04/jimmunol.1700893)

* add path of the program to $PATH environment variable:

  export PATH=$PATH:/path/to/netMHCpan

* download CassOpt using git:

  git clone https://github.com/open-projects/CassOpt

* test CassOpt:

  cd ./CassOpt

  CassOpt.py -f ./test/input_file.fa
  
* you can specify the path to netMHCpan4 program using the option -p:
  
  CassOpt.py -p /path/to/netMHCpan4 -f ./test/input_file.fa
  


## Documentation

Detailed description of CassOpt can be found in the [manual](https://github.com/open-projects/CassOpt/blob/master/UserManual.pdf)

If you haven't found the answer to your question in the docs, or have any suggestions concerning new features, feel free to create an issue here, on GitHub, or write an email to dmitry.malko at gmail.com:
<br />![my mail](https://user-images.githubusercontent.com/5543031/28415000-8bea641e-6d56-11e7-85ca-4287500a4192.png)

## License
Copyright (c) 2018, 2019, D. Malko
All Rights Reserved

PeptoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).



