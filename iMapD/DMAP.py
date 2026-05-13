import numpy as np
import mdtraj as mdt
import tqdm

from scipy.linalg import eigh,eig


def RMSDMatCalc(dataFile = None,
               topFile = ''):
  
  
    if dataFile is None:
        print('Data is not provided')
        KeyboardInterrupt
    if topFile == '':
        print('Topology is not provided')
        KeyboardInterrupt
    
            

    trajectory = mdt.load(dataFile,top=topFile)

    rmsdMatrix = np.zeros((trajectory.n_frames,trajectory.n_frames))

    for i in tqdm.tqdm(range(trajectory.n_frames-1)):
        rmsdMatrix[i,i+1:] = mdt.rmsd(trajectory[i+1:]
                                        ,trajectory[i]
                                        ,atom_indices=trajectory.topology.select('backbone'))*10

    rmsdMatrix = rmsdMatrix+rmsdMatrix.T
    
    return rmsdMatrix
  
  
def dmapFunc(distanceMat=None, epsilon=0.0):
    
    if distanceMat is None :
        print('distanceMat is not provided')
        raise KeyboardInterrupt
    if epsilon == 0.:
        print('epsilon (scaling parameter) is not given')
        raise KeyboardInterrupt
    
    
    KMat = np.exp(-distanceMat**2 / epsilon**2)
    
    
    r = np.sum(KMat, axis=0)
    Di = np.diag(1/r)
    P = np.matmul(Di, KMat)
    
    D_right = np.diag((r)**0.5)
    D_left = np.diag((r)**-0.5)
    P_prime = np.matmul(D_right, np.matmul(P,D_left))
 
    eigenValues, eigenVectors = eigh(P_prime)
    
    
    ## Sorting the components
    sortingIndecies = eigenValues.argsort()[::-1]
    eigenValues = np.real(eigenValues[sortingIndecies])
    eigenVectors = np.real(eigenVectors[:,sortingIndecies])
    
    diffusion_coordinates = np.matmul(D_left, eigenVectors)
    
    return eigenValues,diffusion_coordinates
