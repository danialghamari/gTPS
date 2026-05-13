import numpy as np
import os
import sys

from simulation_utilities import *
from DMAP import *
from NewPoint import *
from analysis import *

import mdtraj as mdt

from scipy.spatial import ConvexHull

import multiprocessing



'''
    The main part of code for running iMapD simulation. The full description is given in : https://doi.org/10.1073/pnas.1621481114
'''
numOfThreads = multiprocessing.cpu_count()
recursiveNumSimulation = 5 ## Number of times to resample and run MD from a new point, generated outside of previously explored region
currentDir = os.getcwd()

topFile = currentDir+'/path/to/initial/config'

os.system('mkdir -p DATA/InitialSimulation')

gmx_ob = Gromacs()

os.chdir('DATA/InitialSimulation')

prep4Sim(gmx_ob=gmx_ob , initconf=topFile , mdpfiles_path='/path/to/mdp')
gmx_ob.gmxcmd(cmd = 'grompp',f = '/path/to/mdp/MD.mdp'
                    ,c = 'NPT.gro'
                    ,r = 'NPT.gro'
                    ,t = 'NPT.cpt'
                    ,p = 'topol.top'
                    ,o = 'NPT.tpr')
gmx_ob.gmxcmd('v',cmd = 'mdrun' ,deffnm = 'MD')

gmx_ob.gmxcmd('center', cmd = 'trjconv', 
                needs_input=True, _input = ('1','1'),
                s = 'MD.tpr',
                f = 'MD',
                ur = 'compact',
                pbc = 'mol',
                o = 'MD_noPBC.xtc')

os.chdir(currentDir)

topFile = currentDir+'/DATA/InitialSimulation/initConf.gro'
nptFile = currentDir+'/DATA/InitialSimulation/MD.gro'
allDataFile = currentDir+'/DATA/alldata.xtc'

data = traj2Numpy(dataPath = allDataFile
                ,topFile= topFile)

numpy2PDB(data = data ,dirSave= allDataFile ,topFile=topFile ,nptFile= nptFile)

for i in range(int(sys.argv[1]),int(sys.argv[2])):
    
    print('ITERATION {:d}'.format(i+1))
    iterDir = currentDir+'/DATA/Iter_{:d}'.format(i+1)

    os.system('mkdir -p '+iterDir)
    
    numpy2PDB(data = data ,dirSave= allDataFile ,topFile=topFile ,nptFile= nptFile)
    
    rmsdMatrix = RMSDMatCalc(dataFile=allDataFile, topFile= topFile)

    minimumDistance = []
    for i,j in enumerate(rmsdMatrix[:-1]):
        minimumDistance.append(np.sort(rmsdMatrix[i])[1])


    eigval, dmapComps = dmapFunc(rmsdMatrix,np.mean(rmsdMatrix[rmsdMatrix>0]))

    dmapComps = dmapComps[:,1:3]

    convexVertices = ConvexHull(dmapComps).vertices
    
    
    newPointsDict,neighborsDict = NewPoints(convexVertices
                                ,data=data
                                ,distanceMat=rmsdMatrix
                                ,const = 1.
                                ,tolDistance= np.max(minimumDistance)
                                ,tolNumNeighbors= 30 # minimum 30 neighbors
                                ,initShape=(data.shape[1],data.shape[2]))


    for j in newPointsDict:
        print('\t POINT {:d}'.format(j+1))
        
        pointDir = iterDir+'/Point_{:d}'.format(j+1)
        
        os.system('mkdir -p '+pointDir+'/Ratchet')
        
        initPointFile = iterDir+'/Point_{:d}'.format(j+1)+'/initialPoint.pdb'
        targetPointFile = iterDir+'/Point_{:d}'.format(j+1)+'/target.pdb'
        
        
        numpy2PDB(data = neighborsDict[j] ,dirSave= initPointFile ,topFile=topFile ,nptFile= nptFile)
        numpy2PDB(data = newPointsDict[j] ,dirSave= targetPointFile ,topFile=topFile ,nptFile= nptFile)

        os.chdir(pointDir+'/Ratchet')
        Structure(targetPointFile).save_sarr('target.cmp')

        prep4Sim(gmx_ob=gmx_ob , initconf=targetPointFile , mdpfiles_path='/path/to/mdp')
        gmx_ob.gmxcmd(cmd = 'grompp',f = '/path/to/mdp/rMD.mdp'
                            ,c = 'NPT.gro'
                            ,r = 'NPT.gro'
                            ,t = 'NPT.cpt'
                            ,p = 'topol.top'
                            ,o = 'NPT.tpr')
        gmx_ob.gmxcmd('v',cmd = 'mdrun' ,deffnm = 'MD')

        gmx_ob.gmxcmd('center', cmd = 'trjconv', 
                        needs_input=True, _input = ('1','1'),
                        s = 'MD.tpr',
                        f = 'MD',
                        ur = 'compact',
                        pbc = 'mol',
                        o = 'MD_noPBC.xtc')

        os.chdir(currentDir)

        
        ratchetData = traj2Numpy(dataPath = pointDir+'/Ratchet/MD_noPBC.xtc'
                ,topFile= topFile)
        
        newPointFile = pointDir+'/point.gro'
        nptFile = pointDir+'/Ratchet/MD.gro'
        
        numpy2PDB(data = ratchetData[-1],
                 dirSave = newPointFile,
                 topFile = topFile,
                 nptFile = nptFile)
        
        del ratchetData
        
        for k in range(recursiveNumSimulation):
            print('\t\t ITER {:d}'.format(j+1))

            newDataDir = iterDir+'/Point_{:d}'.format(j+1)+'/Iter_{:d}'.format(k+1)

            os.mkdir(newDataDir)

            os.chdir(newDataDir)
            prep4Sim(gmx_ob=gmx_ob , initconf=newPointFile , mdpfiles_path='/path/to/mdp')
            gmx_ob.gmxcmd(cmd = 'grompp',f = '/path/to/mdp/consecutive_MD.mdp'
                                ,c = 'NPT.gro'
                                ,r = 'NPT.gro'
                                ,t = 'NPT.cpt'
                                ,p = 'topol.top'
                                ,o = 'NPT.tpr')
            gmx_ob.gmxcmd('v',cmd = 'mdrun' ,deffnm = 'MD')

            gmx_ob.gmxcmd('center', cmd = 'trjconv', 
                            needs_input=True, _input = ('1','1'),
                            s = 'MD.tpr',
                            f = 'MD',
                            ur = 'compact',
                            pbc = 'mol',
                            o = 'MD_noPBC.xtc')

            os.chdir(currentDir)
            
            data = np.concatenate((data,traj2Numpy(dataPath = newDataDir+'/MD_noPBC.xtc'
                                                   ,topFile= topFile)),axis=0)
            
    
    # ### keeping a clean directory for the data
    # os.chdir(iterDir)
    # os.system("find . -type f -not \( -name *_noPBC.xtc -o -name *point* -o -name target.pdb -o -name initial* -o -name initConf.gro -o -name mdout.mdp -o -name ratchet.out -o -name *Ratchet* \) -delete")
    # os.chdir(currentDir)
        
            
numpy2PDB(data = data ,dirSave= allDataFile ,topFile=topFile ,nptFile= nptFile)
