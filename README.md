# gTPS: graph-Transition Path Sampling
The code I developed and/or implemented during my PhD research. The related articles: [DOI:10.1038/s41598-022-20032-x](https://doi.org/10.1038/s41598-022-20032-x) and [DOI:10.1021/acs.jctc.3c01174](https://doi.org/10.1021/acs.jctc.3c01174)

The main Python codes are the implementation of [iMapD framework](#imapd) which performs a manifold search (for macromolecular systems) through usage of [DiffusionMap](#dmap) to parameterize the explored parts of underlying intrinsic manifold and then apply [PCA](#pca) (Principal Component Analysis) to generate a data point outside. The frameowrk had shown up to 100x speed-up over plain molecular dynamics.

The other files are Jupyter Notebooks which we used to perform clustering, building the network of connectivity, and eventually obtain pathways using Monte Carlo and a Quantum Annealer (D-Wave).

The external packages:
Numpy\
Scipy\
Scikit-learn\
Matplotlib\
MDAnalysis\
MDTraj\
OceanSDK (DWave)


<a name="imapd"></a>[DOI:10.1073/pnas.1621481114](https://doi.org/10.1073/pnas.1621481114)\
<a name="dmap"></a>[DOI:10.1016/j.acha.2006.04.006](https://doi.org/10.1016/j.acha.2006.04.006)\
<a name="pca"></a>[DOI:10.1145/3447755](https://doi.org/10.1145/3447755)
