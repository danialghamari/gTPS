import numpy as np
from scipy.spatial import ConvexHull

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as skPCA


##Finding neighboring clusters based on Euclidean distance, for boundary points
def NeighborIdentifier(boundaryPointsIndecies = None
                      , data = None
                      , distanceMat = None
                      , tolDistance = 0.
                      , tolNumNeighbors=0.):
    
            clustersDict = {}

            for i,j in enumerate(boundaryPointsIndecies):

                distanceArray = distanceMat[j]
                ## we do not include the boundary points here
                indexBoundaryNeighbors = np.where((distanceArray <= tolDistance) &\
                                                 (distanceArray > 0))[0]

                ## if each cluster has less than tolNeighbors points, we increase it to tolNeighbors
                if indexBoundaryNeighbors.shape[0] < tolNumNeighbors:
                    indexBoundaryNeighbors = np.argsort(distanceArray)[1:tolNumNeighbors]
                
                ## we put the boundary point in the cluster
                clustersDict.update( { \
                                    i:np.concatenate((data[j].reshape(1,-1),\
                                                      data[indexBoundaryNeighbors].reshape(indexBoundaryNeighbors.shape[0],-1)),\
                                                      axis=0 ) } )

            return clustersDict


### apply pca on cluster + extension
def NewPoints(boundaryPointsIndecies = None
               , data = None
               , distanceMat = None
               , const = 0. ## distance to extend beyond boundary
               , tolDistance = 0. ## proximity tolerance for each boundary point
               , tolNumNeighbors = 0  ## minimum number of neighbors for each boundary points 
                                      ## (it also determines PCA n_features, should be greater or equal to "dimension of data"
                                      ## if you don't include boundary in the evaluation of PCA, otherwise "dimension of data-1")
                                      ##
            #    , tolPCAComponents = 0 ## number of PCA components required
               , initShape = None  ): ## dimension of data
    
            
    if boundaryPointsIndecies is None:
        print('boundaryPointsIndecies is not given')
        raise KeyboardInterrupt
    if data is None:
        print('data is not given')
        raise KeyboardInterrupt
    if distanceMat is None:
        print('distanceMat is not given')
        raise KeyboardInterrupt
    if const == 0.:
        print('constant is zero')
        raise KeyboardInterrupt
    if tolDistance == 0.:
        print('distance tolerance is zero')
        raise KeyboardInterrupt
    if initShape is None:
        print('initShape is not defined')
        raise KeyboardInterrupt


    clustersNeighborsDict = NeighborIdentifier( boundaryPointsIndecies,\
                                                data,\
                                                distanceMat,\
                                                tolDistance,\
                                                tolNumNeighbors
                                               )

    newPointsDict = {}

    PCA = skPCA()
    SS = StandardScaler()

    for i in clustersNeighborsDict:

        # if clustersNeighborsDict[i].shape[0]<tolNumNeighbors:
        #     Exception('number of neighbors is less than tolNumNeighbors')

        SS.fit(clustersNeighborsDict[i])

        standardizedData = SS.transform(clustersNeighborsDict[i][1:])

        pca = PCA.fit(standardizedData)

        for j in range(pca.explained_variance_ratio_.shape[0]):
            if np.sum(pca.explained_variance_ratio_[:j+1]) > 0.98:
                nComponents = j+1
                break

        loadingMat = pca.components_[:nComponents].T


        yB = clustersNeighborsDict[i][0].dot(loadingMat)

        xCOM = np.mean(clustersNeighborsDict[i][1:],axis=0)

        yCOM = xCOM.dot(loadingMat)

        yNew = yB-yCOM \
                + const*(yB-yCOM)/np.linalg.norm(yB-yCOM)

        newPointsDict.update({ i:(yNew.dot(loadingMat.T) \
                       + xCOM).reshape(initShape) })

    
    neighborsDict = {}
    for i in clustersNeighborsDict:
        neighborsDict.update({i:clustersNeighborsDict[i][0].reshape(initShape)})
    
    return newPointsDict,neighborsDict