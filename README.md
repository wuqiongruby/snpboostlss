---
editor_options: 
  markdown: 
    wrap: 72
---

# snpboostlss

An algorithm to apply statistical boosting for distributional regression
on genotype data via a batch-based approach.

snpboostlss includes R functions in which the boosting algorithm is
implemented.

Simulations provides R Code and bash scripts that were used to run the
simulation on phenotypes based on HAPNEST synthetic genotype data and
fit statistical boosting models on those.

UK_biobank_application provides R scripts that were used to fit
snpboostlss and snpboost models on UKBB data.

Installation: The following requirements of snpboostlss are available
from CRAN:

-   tidyverse

-   data.table

-   Rfast

-   parallel

-   dataPreparation

-   bigsnpr

-   cowplot

-   readr

-   gamlss

-   pgenlibr

Notice that if you are using Windows system, installation of sed is
necessary. sed can be downloaded from
<https://gnuwin32.sourceforge.net/packages/sed.htm>.

Additionally PLINK 2.0 is required. It can be installed from
<https://www.cog-genomics.org/plink/2.0/>.
