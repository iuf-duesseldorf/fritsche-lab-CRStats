# CRStats - Concentration-Response Statistics

A R package for concentration-response analysis automation for in vitro test systems optimized for multi-well plate experiments based on drc.

### Package Ecosystem

**CRStats** builds upon the following packages not developed by members of the **CRStats** team:

- R6: Reference class objects.
- data.table: Extension of R's data.frame.
- drc: Analysis of Dose-Response Curves

All these packages are well curated and mature; we expect no problems with dependencies.
Additionally, we suggest the following packages for extra functionality:

- For parallelization: parallel.
- For progress bars: progressr.
- For capturing output, warnings, and exceptions: evaluate or callr.


## Team
### Developer

**Arif Dönmez**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

### Correspondence

**Ellen Fritsche**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

**Martin Scholze**: College of Health, Medicine and Life Sciences\
Brunel University London, Uxbridge, UK

**Eike Keßel**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

**Stefan Masjosthusmann**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

**Katharina Koch**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

**Kristina Bartmann**: Alternative method development for environmental toxicity testing\
IUF - Leibniz Research Institute for Environmental Medicine

## Installation
### Alpha version
``` r
## You can install CRStats from GitHub
# install.packages("devtools")
#
devtools::install_github("iuf-duesseldorf/fritsche-lab-CRStats")
```
### Latest developer editon with experimental features
``` r
## You can install CRStats from GitHub
# install.packages("devtools")
#
devtools::install_github("ArifDoenmez/CRStats")
```
