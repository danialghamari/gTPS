import mdtraj as mdt
import numpy as np
from numba import njit, float32, prange


@njit("float64[:,:](float32[:,:])")
def CArrayFunc(distances = None):
    
    c_array = np.zeros(distances.shape)
    
    beta = 5.0
    constant = 1.5
    r0 = 0.75
    
    c_array = 1/(1+np.exp(beta*(distances-constant*r0)))

    return c_array


@njit("float64[:,:](float64[:,:])",parallel=True)
def ContactDistanceFunc(c_array):
    
    contactDistance_matrix = np.zeros((c_array.shape[0],c_array.shape[0]))
    
    for i in prange(c_array.shape[0]-1):
        contactDistance_matrix[i,i+1:] = np.sqrt(np.sum(np.power(c_array[i+1:]-c_array[i],2),axis=1))
        
    
    contactDistance_matrix = contactDistance_matrix+contactDistance_matrix.T

    return contactDistance_matrix



def contactDistances(trajectory = None):

    if trajectory is None:
        print('trajectory must be defined')
        raise KeyboardInterrupt

    
    distances, res_pairs = mdt.compute_contacts( trajectory, contacts = "all"
                                                , scheme = "ca"
                                                , ignore_nonprotein = True
                                                , periodic = False )

    
    TRAJ_c_array = CArrayFunc(distances)

    
    TRAJ_contactDistance_matrix = ContactDistanceFunc(TRAJ_c_array)
    
    return TRAJ_contactDistance_matrix