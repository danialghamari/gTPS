# gTPS: graph-Transition Path Sampling
The code I developed and/or implemented during my PhD research. The related articles: \href{doi.org/10.1038/s41598-022-20032-x}{DOI:10.1038/s41598-022-20032-x} and \href{https://doi.org/10.1021/acs.jctc.3c01174}{DOI:10.1021/acs.jctc.3c01174}

The main Python codes are the implementation of iMapD framework[^imapd] which performs a manifold search (for macromolecular systems) through usage of DiffusionMap[^dmap] to parameterize the explored parts of underlying intrinsic manifold and then apply PCA[^PCA] (Principal Component Analysis) to generate a data point outside. The frameowrk had shown up to 100x speed-up over plain molecular dynamics.

The other files are Jupyter Notebooks which we used to perform clustering, building the network of connectivity, and eventually obtain pathways using Monte Carlo and a Quantum Annealer (D-Wave).

The external packages:
Numpy
Scipy
Scikit-learn
Matplotlib
MDAnalysis
MDTraj
OceanSDK (DWave)

[^imapd] \href{https://doi.org/10.1073/pnas.1621481114}{DOI:10.1073/pnas.1621481114}
[^dmap] \href{https://doi.org/10.1016/j.acha.2006.04.006}{DOI:10.1016/j.acha.2006.04.006}
[^PCA] \href{https://doi.org/10.1145/3447755}{DOI:10.1145/3447755}
