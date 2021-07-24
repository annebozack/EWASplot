# EWASplot

Common functions for plotting results of epigenome-wide association studies (EWAS) implemented in R and based on ggplot2(https://ggplot2.tidyverse.org/). 

To install:
    `devtools::install_github("annebozack/EWASplot")`

### Q-Q plot example
```
qq_plot(probes$P.Value)
```

![qq plot](https://raw.githubusercontent.com/annebozack/images/master/qq_ex_sm.png)
