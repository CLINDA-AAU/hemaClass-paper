# [hemaClass.org](http://hemaclass.org)
This is the reprository for the paper presenting the **hemaclass** [package](https://github.com/oncoclass/hemaclass) and [website](http://hemaclass.org)  for classification of diffuse large B-cell lymphomas (DLBCL).
It proposes general microarray preprocessing techniques to, oxymoronically, perform one-by-one robust multichip average (RMA) preprocesssing opposed to regular RMA which needs an entire cohort of patients or samples.
This effort is motivated by personalized medicine as well as other scientific deciplines in which cohort-based RMA normalization is infeasible.
The package and website implements various classifiers of DLBCL [references].

## How to recreate the paper
Clone this github repository and compile the `hemaClass.tex` file with `pdflatex`using the `Makefile` supplied. From the command line or terminal, simply write 
```sh
make
```
to build the paper.

To run the analysis, generating the figures and tables, and build the paper execute
```sh
make all
```
from the command line or terminal. **Warning**, the R-script `hemaclass.R` automatically downloads all the microarray data and performs various data preprocessing. The data needs about 42 GBs of free disk space and the preprocessing is quite memory intensive.

## How to cite this paper
Citation information goes here

## Usage

## References
[References goes here. BAGS and other relevant references.]
