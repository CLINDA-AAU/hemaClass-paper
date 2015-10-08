# [hemaClass.org](http://hemaclass.org)
This is the reprository for the paper presenting the **hemaclass** [package](https://github.com/oncoclass/hemaclass) and [website](http://hemaclass.org)  for classification of diffuse large B-cell lymphomas (DLBCL).
It proposes general microarray preprocessing techniques to, oxymoronically, perform one-by-one robust multichip average (RMA) preprocesssing opposed to regular RMA which needs an entire cohort of patients or samples.
This effort is motivated by personalized medicine as well as other scientific deciplines in which cohort-based RMA normalization is infeasible.
The package and website implements various classifiers of DLBCL [1].


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


## Usage
A tutorial of how to use the website and/or package is provided on the [website](http://hemaclass.org) and as a vignette in the package, e.g. run `vignette(howto, package = "hemaClass")`.


## Citing this paper
We have invested alot of time and effort into the paper and [**hemaClass**](https://github.com/oncoclass/hemaclass), so please cite them if you use elements of them. Please confer `citation("hemaClass")`.

1. Steffen Falgreen, Anders Ellern Bilgrau, Jonas Have; **"hemaClass: Online classification of gene expression profiles in hematological cancers."** (2015) R package version 0.11; http://github.com/falgreen/hemaClass

2. Steffen Falgreen, Anders Ellern Bilgrau, Lasse Hjort Jakobsen, Jonas Have, Kasper Lindblad Nielsen, Tarec Christoffer El-Galaly, Julie Støve Bøker, Alexander Schmitz, Hans Erik Johnsen, Karen Dybkær, and Martin Bøgsted; **"hemaClass.org: An online based diffuse large B-cell lymphoma classification tool."** (2015) In preperation for BMC genomics.



## References

1. Dybkær K, Bøgsted M, Falgreen S, Bødker JS et al. *"Diffuse Large B-cell Lymphoma Classification System That Associates  Normal B-cell Subset Phenotypes with Prognosis."* Journal of Clinical Oncology 33, no. 12 (2015): 1379-1388. (GSE56315)
       
2. Falgreen S, Dybkær K, Young KH, Xu-Monette ZY et al. *"Predicting response to multidrug regimens in cancer patients using cell line experiments and regularised regression models."* BMC cancer 15, no. 1 (2015): 235.

3. Laursen, MB, Falgreen S, Bødker JS, Schmitz A, et al. *"Human B-cell cancer cell lines as a preclinical model for studies of drug effect in diffuse large B-cell lymphoma and multiple myeloma."* Experimental Hematology 42, no. 11 (2014): 927-938.
