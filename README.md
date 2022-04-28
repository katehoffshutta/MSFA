# MSFA

author: Roberta De Vito, Ruggero Bellio

this fork is modified by: Kate Hoff Shutta (kshutta@hsph.harvard.edu)

Fits the Multi-Study Factor Analysis model  via  [ECM algorithm](#1-fitting-a-msfa-model-via-the-ecm-algorithm) and perform [graphical modeling](#2-graphical-modeling) of multi-study data with an extension.

## 1 Fitting a MSFA model via the ECM Algorithm

The following example illustrates how to fit a MSFA model via the ECM Algorithm, 
using a data set available in the Bioconductor repository (www.bioconductor.org). 

### Getting the data
Some pre-processing is required to get the data into a form suitable for the
analysis. This was already done, and the resulting data frame is saved into the
`data_immune` object. The commands that were used to form it are included
in the help file for the data object.

```{r help2, echo = TRUE, results = TRUE, tidy = TRUE}
library(MSFA)
data(data_immune)
help(data_immune)
```

### Obtaining suitable starting values for model parameters
Then we get suitable starting values for model parameters, selecting K=3 common
factors and (3, 4) study-specific factors for the two studies, respectively.

```{r, starting values, messages = FALSE}
start_value <- start_msfa(X_s = data_immune, k = 3, j_s = c(3, 4))
```

### Fitting the model via ECM
Now everything is in place for estimating the model parameters via the ECM algorithm

```{r get estimate, results = FALSE}
mle <-  ecm_msfa(data_immune, start_value, trace = FALSE)
```

The estimated matrix of common loadings can be represented by a suitable heatmap:

```{r heatmap, fig.show = 'hold', fig.width = 7.5, fig.height = 6.5, message = FALSE}
library(gplots)
heatmap.2(mle$Phi,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', density.info="none", col=heat.colors(256))
```

---
 ## 2 Graphical Modeling

Under construction.

