# OPWpaper1
An R package to reproduce the results and figures from simulations and data described in the Optimal Pvlaue weighting Paper

# Installing the package and reproducing all simulations and data examples

You can install the package as follows:

```{r}
library("devtools")
# install OPWeight
install_github("vladchimescu/lpsymphony", subdir="lpsymphony")
install_github("mshasan/OPWeight")

# Bioconductor prerequisites
source("http://bioconductor.org/biocLite.R")
biocLite(c("genefilter","DESeq2","qvalue","Biobase",
            "BiocParallel","airway","pasilla", "BiocStyle"))
# finally install this package
install_github("mshasan/OPWpaper1")
```

