# EWASplot

Common functions for plotting results of epigenome-wide association studies (EWAS) implemented in R and based on ggplot2(https://ggplot2.tidyverse.org/). 

To install:
    `devtools::install_github("annebozack/EWASplot")`

### Q-Q plot 

Quantile-quantile plot of observed and expected p-values. Genomic inflation factor (λ) is calculated using the median method.

```
qq_plot(probes$P.Value)
```

![qq plot](https://raw.githubusercontent.com/annebozack/images/master/qq_ex_sm.png)

### Volcano plot 

Volcano plot of effect estimates and p-values. Aesthetics including point size and color can be modified.

```
volcano_plot(probes$P.Value, probes$logFC, FDR = T)
```

![vol plot](https://raw.githubusercontent.com/annebozack/images/master/vol_ex.png)

### Manhattan plot 

Manhattan plot of p-values and genomic locations. If results from 450K or EPIC analyses are provided with genomic locations, probes will be annotated using Bioconductor annotation packages. The locations of regions can also be plotted by providing a region argument. Aesthetics including point size and color can be modified.

```
manhattan_plot(probes, regions, FDR = T)
```

![manhattan plot](https://raw.githubusercontent.com/annebozack/images/master/man_ex_sm.png)

