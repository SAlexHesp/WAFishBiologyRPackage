# WAFishBiologyRPackage

Collection of methods for analysis of fish biological data.
 
"WAFishBiology" contains a range of analyses associated with estimating important biological parameters typically required for stock assessments from fish/invertebrate data, including the parameters of various growth, maturity and weight-length parameter relationships. The package also contains some other analyses and plotting functions for exploring fish reproductive data, such as for monthly gonadosomatic indices and monthly proportions of fish with gonads at different developmental stages. The analyses in this package are also intended to complement the R package "L3Assess", which contains a range of age and length-based catch curve and per recruit anlayses, applicable for data-limited assessments. 

As with other packages in github, to install WAFishBiology, first ensure you have the 'devtools' package, otherwise install using install.packages("devtools').
Then, use: 

library(devtools)

devtools::install_github("SAlexHesp/WAFishBiologyRPackage", build_vignettes=TRUE)

If you already have a version of the WAFishBiology package installed but wish to update, I suggest using the line of code below rather than the one above 
to ensure the updated version is installed.
devtools::install_github("SAlexHesp/WAFishBiologyRPackage", build_vignettes=TRUE, force=TRUE)
