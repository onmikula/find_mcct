### TOOLS
library(ape)
library(phangorn)
source("func/find_mcct.R")


### ultrametric trees from species tree inference in StarBEAST2 (Ogilvie et al. 2017)
# species trees of African pygmy mice, subgenus Nannomys (Krásová et al. 2022)
nann <- ape::read.tree("data/Nannomys.trees")
# maximum clade credibility tree with mean common ancestor node heights
mcct_comanc <- find_mcct(trees=nann, file="results/Nannomys_mcct.tree", monophyletic=FALSE, method="mean", return=TRUE)
# maximum clade credibility tree with median node heights calculated from a subset of trees where the clades are monophyletic (Drummond & Bouckaert 2015, pp.)
mcct_median <- find_mcct(trees=nann, file="results/Nannomys_mcct.tree", monophyletic=TRUE, method="median", return=TRUE)
# the second most probably topology (M. "Nyika" sister to the rest) used as a target to calculate mean common ancestor node heights
nyika <- "((((((tritonA,tritonB),tritonCD),Harenna),imberbis),((neavei,(((gerbillus,minutoides),gratus),mattheyi)),(bufo,(proconodon,mahomet)))),Nyika);"
target_comanc <- find_mcct(trees=nann, file="results/Nannomys_nyika.tree", target=nyika, monophyletic=FALSE, method="mean", return=TRUE)


### unrooted trees inferred in MrBayes (Ronquist et al. 2012)
# cytochrome b tree of M. triton and M. harennensis with three non-African Mus as outgroups (Krásová et al. 2022)
# posterior samples from two independent runs
cytb1 <- ape::read.nexus("data/Mharennensis_cytb.run1.t")
cytb2 <- ape::read.nexus("data/Mharennensis_cytb.run2.t")
# discard burnins (here specified as percentage) and combine the posterior samples
cytb <- combine_runs(list(cytb1, cytb2), burnin=10)
# maximum bipartition credibility tree (unrooted), argument monophyletic ignored as monophyletic=FALSE is not defined
mbct_mean <- find_mcct(trees=cytb, file="results/Mharennensis_cytb_mbct.tree", method="mean", return=TRUE, thinning=5)
# use outgroups to root the trees and drop them
outgroups <- c("AB125782_Mplatythrix", "AB819919_Mmusculus", "AB125769_Mcookii")
cytb <- root_and_drop(cytb, outgroup=outgroups)
# maximum clade credibility tree with mean common ancestor node heights calculated from rooted, but non-ultrametric trees
mcct_nonultr <- find_mcct(trees=cytb, file="results/Mharennensis_cytb_mcct.tree", monophyletic=FALSE, method="mean", return=TRUE)


### time-calibrated but non-ultrametric trees from fossilized birth-death dating in BEAST2 (Gavryushkina et al. 2017)
# concatenated trees of Bathyergidae (African mole-rats) and related phyomyid rodents, including dated fossils (Uhrová et al. 2022)
bathy <- ape::read.nexus("data/Bathyergidae_dating.trees")
# fossil species removed
bathy <- remove_fossils(bathy, tol=1e-5)
# maximum clade credibility tree
# 20% burnin discarded (could be done also prior to mcc tree estimation by 'discard_burnin' function)
# every 5th tree taken for the calculation
mcct_nofossils <- find_mcct(trees=bathy, file="results/Bathyergidae_dating_mcct.tree", monophyletic=FALSE, method="mean", return=TRUE, burnin=0.20, thinning=5)


### REFERENCES
#Gavryushkina A, Heath TA, Ksepka DT, Stadler T, Welch D, Drummond AJ (2017) Bayesian total-evidence dating reveals the recent crown radiation of penguins. Systematic Biology, 66: 57-73.  https://doi.org/10.1093/sysbio/syw060.
#Krásová J, Mikula O, Lavrenchenko LA, Šumbera R, Meheretu Y, Bryja J (2022) A new rodent species of the genus Mus (Rodentia: Muridae) confirms the biogeographical uniqueness of the isolated forests of southern Ethiopia. Organisms Diversity & Evolution, accepted.
#Ogilvie HA, Bouckaert RR, Drummond AJ (2017) StarBEAST2 brings faster species tree inference and accurate estimates of substitution rates. Molecular Biology & Evolution, 34: 2101-2114.  https://doi.org/10.1093/molbev/msx126
#Ronquist F, Teslenko M, Van Der Mark P, Ayres DL, Darling A, Höhna S, Larget B, Liu L, Suchard MA, Huelsenbeck JP (2012) MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic Biology, 61: 539-542. https://doi.org/10.1093/sysbio/sys029
#Uhrová M, Mikula O, Bennett NC, Van Daele P, Piálek L, Bryja J, Visser JH, van Vuuren BJ, Šumbera R (2022) Species limits and phylogeographic structure in two genera of solitary African mole-rats Georychus and Heliophobius. Molecular Phylogenetics and Evolution, 167: 107337. https://doi.org/10.1016/j.ympev.2021.107337



