# find_mcct
R function post-processing outputs of phylogenetic placement analysis.

**find_mcct** is an R function (R Core Team 2021) which employs 'maxCladeCred' of 'phangorn' (Schliep 2011) to find maximum clade credibility tree with specified type of node heights (and hence also branch lengths). The function also depends on R package ape (Paradis & Schliep 2019).

For more details about the function and its arguments, see comments in the file func/find_mcct.R and in the example routines in find_mcct_example.R. 

**References:**
- Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526–528. https://doi.org/10.1093/bioinformatics/bty633
- R Core Team (2021) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
- Schliep KP (2011) phangorn: phylogenetic analysis in R. Bioinformatics, 27: 592-593. https://doi.org/10.1093/bioinformatics/btq706
