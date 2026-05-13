
from analysis import *
import MDAnalysis as mda, os, mdtraj as mdt

'''
    An API class to run GROMACS commands in python scripts
    Similar to running GMX from command-line.
    
    Boolean flags must be given as *args and before any other non-boolean flags
    Non-boolen flags, e.g -f /some/file.gro, must be given as kwargs, e.g f="/some/file"

    Depending on the HPC configuration you might want to redefine how GROMACS excutes the "mdrun" command (MPI or OpenMP)
    Here, it runs with MPI config
'''


class Gromacs():
    def __init__(self,
                 np_mpiexec= 1, ## 
                 ntomp=4):
        self.np_mpiexec =  np_mpiexec
        self.ntomp = ntomp

    def set_mpi_omp(self,
                    mpi_=2,
                    omp_=1):
        self.np_mpiexec = mpi_
        self.ntomp = omp_

    def gmxcmd(self, 
               *args,
               cmd = '', ## The command to be executed
               needs_input = False, ## In case the command requires an input. The input will be passed as : printf 'input1\ninput2\n...' | gmx_mpi ...
               _input = (), ## If needs_input = True, you need to provide the inputs 
               **kwargs):
        if cmd=='':
            raise ValueError("please provide a command")
        if needs_input and len(_input)==0:
            raise ValueError("please provide the input for the command, e.g for genion one needs SOL")
            
        
        locals_args = locals()
        
        passing2GMX_args = {}
        passing2GMX_input = 'printf \''+''.join([f'{inp}{"\\n"}' for inp in _input])+'\'| '
        
        if 'kwargs' in  locals_args.keys():
            passing2GMX_args = passing2GMX_args | {key:item for key,item in locals_args['kwargs'].items()}
        if 'args' in  locals_args.keys():
            passing2GMX_args = passing2GMX_args | {item:'' for item in locals_args['args']}
            
        if cmd!='mdrun':
            os.system(passing2GMX_input if needs_input else '' + f'gmx_mpi {cmd} '+''.join([f' -{key} {item} ' for key,item in passing2GMX_args.items()]))
        else:
            os.system(f'mpiexec -np {self.np_mpiexec} gmx_mpi {cmd} -ntomp {self.ntomp} '+''.join([f' -{key} {item} ' for key,item in passing2GMX_args.items()]))


'''
    A generic function to prepare the system for simulation
    It could be used a template to redefine some steps if needed.
'''
def prep4Sim(gmx_ob, initconf = '', mdpfiles_path=''):

    if initconf == '' or mdpfiles_path == '':
        raise ValueError("please provide an initial config. and path to mdpfiles")
    
    gmx_ob.gmxcmd('ignh',cmd= 'pdb2gmx', f='initConf.gro', o = 'init.gro', ff='amber99sb-ildn', water='tip3p' )
    gmx_ob.gmxcmd(cmd= 'editconf', f = 'init.gro', o='Box.gro', c='yes',d = '1', bt='cubic' )
    gmx_ob.gmxcmd(cmd= 'solvate', cp = 'Box.gro' ,o = 'Solvated.gro' ,cs = 'spc216.gro' ,p = 'topol.top')
    gmx_ob.gmxcmd(cmd= 'grompp', f = mdpfiles_path+'/ions.mdp' ,c = 'Solvated.gro' ,p = 'topol.top' ,o = 'Solvated.tpr')
    gmx_ob.gmxcmd('neutral', cmd= 'genion', s = 'Solvated.tpr' ,o = 'Neutralized.gro' ,p = 'topol.top')

    gmx_ob.gmxcmd(cmd='grompp'
                ,f = mdpfiles_path+'/EnergyMinimization.mdp' 
                ,c = 'Neutralized.gro' 
                ,p = 'topol.top' 
                ,o = 'EnergyMinimize.tpr')
    gmx_ob.gmxcmd('v', cmd='mdrun',deffnm = 'EnergyMinimize')

    gmx_ob.gmxcmd(cmd = 'grompp',f = mdpfiles_path+'/NVT.mdp'
                    ,c = 'EnergyMinimize.gro'
                    ,r = 'EnergyMinimize.gro'
                    ,p = 'topol.top'
                    ,o = 'NVT.tpr')
    gmx_ob.gmxcmd('v',cmd = 'mdrun',deffnm = 'NVT')



    gmx_ob.gmxcmd(cmd = 'grompp',f = mdpfiles_path+'/NPT.mdp'
                    ,c = 'NVT.gro'
                    ,r = 'NVT.gro'
                    ,t = 'NVT.cpt'
                    ,p = 'topol.top'
                    ,o = 'NPT.tpr')
    gmx_ob.gmxcmd('v',cmd = 'mdrun' ,deffnm = 'NPT')



### Importing the data as a Numpy array, the unit is in angstrom
def traj2Numpy(dataPath = ''
                , topFile = ''):
    
    if dataPath == '':
        raise ValueError('Please provide the path to trajectory file')
    if topFile == '':
        raise ValueError('Please provide the topology file')

    mainTraj = mdt.load(dataPath,top=topFile)

    return mainTraj.xyz*10 ### To make it angstrom


### Change last frame of Data to PDB
def numpy2PDB(data = None
            ,dirSave = ''
            ,topFile = ''
            ,nptFile = ''):

    if data is None:
        raise ValueError('Please provide the data')
    if dirSave == '':
        raise ValueError('Please provide the path to save')
    if topFile == '':
        raise ValueError('Please provide the topology file')
    if nptFile == '':
        raise ValueError('Please provide the NPT file')
        
    nptUniverse = mda.Universe(nptFile)
    boxLength = nptUniverse.trajectory[0].dimensions[0]

    universe = mda.Universe(topFile,data)
    
    universe.dimensions = [boxLength,boxLength,boxLength,90,90,90]

    protein = universe.select_atoms('protein')
    
    with mda.Writer(dirSave, protein.n_atoms) as W:
        for ts in universe.trajectory:
            W.write(protein)


