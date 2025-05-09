# ======================================================================
#     Import modules
# ======================================================================
import os, sys, copy, time
# from petsc4py import PETSc
from mpi4py import MPI
import numpy as np
# from mdo_regression_helper import *
from adflowtransi import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem, AeroTransiProblem, TransiProblem
from transi import TransitionCalc
from pyaerotransi import AeroTransi
from multipoint import *
import argparse
from multipoint import createGroups
from pyspline import *
from pygeo import *
# from pywarp import *
# ================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--mach", type=str, default=0.19)
# parser.add_argument("--mach", type=str, default=2.0)
parser.add_argument("--procs", help="number of processors", type=int, default=16)
parser.add_argument('--alpha', type=str, default=0.7498)
parser.add_argument("--reynolds", type=float, default=5.6e6)
args = parser.parse_args()
# ======================================================================
#     Specify parameters for caculation
# ======================================================================
# rst parameters (beg)
# name = 'rae2822' # Tested!
name = 'test'
# grid_file = './603HLFC.cgns'
# grid_file = './rae2822.cgns'
grid_file = 'test_rans.cgns'
# grid_file = './naca0012_rans-L2.cgns'
output_aero = './output_aero'
output_transi = './output_transi'
# angle of attack
alpha = args.alpha
# mach number
mach = args.mach
# Turbulence intensity
Tu = 0.07 # 0.07%
NCF = 5
NTS = 6
# reynolds
reynolds = args.reynolds
# reference area
areaRef = 1.0
# reference chord
chordRef = 1.0
reynoldsLength = 1.0
# temperature K
T = 288
spanDirection = 'z'
# rst parameters (end)
# ======================================================================
#     Create multipoint communication object
# ======================================================================
gcomm = MPI.COMM_WORLD
npTransi = 1
npAero = gcomm.size - npTransi
if(gcomm.rank == 0):
  if(not(os.path.exists('output_aero'))):
    os.mkdir('output_aero')
  if(not(os.path.exists('output_transi'))):
    os.mkdir('output_transi')
nGroup = 1
nProcPerGroup = args.procs
npTransi = 1
npAero = MPI.COMM_WORLD.size - npTransi
# Creat aero/transition comms
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet('cruise', nMembers=nGroup, memberSizes=nProcPerGroup)
gcomm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()
print('gcomm.rank:',gcomm.rank)
# Creat aero/transition comms
comm, flags = createGroups([npTransi, npAero], comm=gcomm)
aeroID = 1
transiID = 0
print('comm:', comm.rank, flags)
# ======================================================================
#     Options Set-up
# ======================================================================
#rst ADflow aero options(beg)
aeroOptions = {
  # I/O Parameters
  'gridFile':grid_file,
  'outputDirectory':output_aero,
  'monitorvariables':['resrho','cl','cpu','cd','cmz','resturb','cdp','cdv'],
  'volumevariables':['resrho','Intermittency', 'cp', 'mach', 'temp', 'rhoe'],
  'surfacevariables':['cp','vx', 'vy','vz', 'mach','cfx', 'mach', 'rho', 'p', 'temp', 'cf', 'yplus','blank'],
  'writeTecplotSurfaceSolution':True,
	# 'solRestart' : True,
	#'restartFile' : './fc_000_vol.cgns',
  # Physics Parameters
  'equationType':'RANS',
  # Solver Parameters
  'CFL':1.5,
  'CFLCoarse':1.25,
  'MGCycle':'SG',
  'liftIndex': 3,
  # 'rkreset': True,
  # 'nrkreset' : 100,
  # ANK Solver Parameters
  'useANKSolver':True,
  'nsubiterturb':25,
  'anksecondordswitchtol':1e-3, # increased for 30 deg
  # 'ankcflfactor': 4.0,
  # 'ankcoupledswitchtol':1e-5,
  # ANK yayun's setting
	'useANKSolver':True,
  # 'ankuseturbdadi':False,
  # 'ankturbkspdebug':True,
	'ankstepfactor' : 0.5,
	#'ankcoupledswitchtol' : 1e-6,
	'ankmaxiter' : 60,
  # NK Solver Parameters
  'useNKSolver':True,
  'nkswitchtol':1e-8,
  'nkadpc':True,
  # 'nkasmoverlap': 3, # for highly parallel
  'nkinnerpreconits': 2,
  'nkjacobianlag': 3,
  'nkouterpreconits': 3,
  'nkpcilufill': 2,
  'nksubspacesize': 100,
  # Termination Criteria
  'L2Convergence':1e-8,
  'L2ConvergenceCoarse':1e-2,
  'nCycles':100,
  'useblockettes':False,
  #Turb stuff
  # 'useqcr':True, # go 2 NASA tmr for more info (closer to exp??)
  # Following 3 give defalt SA-noft2 (fully turb sim, in lit)
  # 'eddyvisinfratio':.210438,
  'useft2SA' : False,
  # 'turbulenceproduction' : 'vorticity',
  # if use transition
  'ntransition':True,
  'transi2dim':True,
  'useintermittency':True, #Use as is 
  'useintermittencygrid':True, # Keep as False...it is better
	'BoundaryLayerThickness':0.05, # 5-10% of chordref
  'usexyzlstate':False, # Keep as is.
	'nbody':False, # Keep as is.
	'nwingtip':False, # Keep as is.
  'transislicesnum':[2],
  # 'rkreset':True,
  # 'nrkreset':100,
  # Adjoint Parameters
  'setMonitor':False,
  'applyadjointpcsubspacesize':15,
  'adjointL2Convergence':1e-8,
  'ADPC':True,
  'adjointMaxIter': 1000,
  'adjointSubspaceSize':150,
  'ILUFill':2,
  'ASMOverlap':1,
  'outerPreconIts':3,
}
#rst ADflow aero options(end)
#rst transi options(beg)
transiOptions = {
  #Relaxation Parameter
  'outputDirectory':output_transi,
  'RelaxTr': 1.0, #0.8 
  'RelaxLen': 1.0, #1.0
  'TrLimit':1.0,  # 1.0
  'isCompressible':False,#Incompressible or Compressible LST icomp = 1 #Can use false for transonic as well !
  'minerror':1.0, # 1.0 default
  'usexyzlstate':False,# Keep as is.
  'transiiteration':False,# Keep as is.
  'sweeple':[0.0],# Keep as is.
  'sweepte':[0.0],# Keep as is.2D case
  'isLaminarTransition':False, # Meaning -> Docs
  'databasetype':'drela',
  'laminarsepType':'cpcrit',
  'upperfixed':[0.0],
  'lowerfixed':[0.0],
  # 'isquasi3d':False,
  # 'c1factor':0.9,
  # 'usec1avg':True,
  'eNType':0 #0 means LST, 1 mneans database ... keep as 0
  }
#rst transi options(end)
#rst aerotransi options(begin)
lfOptions ={
  'adjointRelTol':1e-7,
  'outputDir':output_aero,
  'adjointDamp0':1.0,
  'adjointsolver':'KSP',
  'lfsubspacesize':40,
  'nadjointiter':100,
  'reltol':1e-3,
  'printtiming':True,
  'nlfiter':3,
  'damp0':1.0
}
#rst aerotransi options(end)
#rst mesh options(begin)
meshOptions = {
  'gridFile': grid_file,
  }
#rst mesh options(end)
# ======================================================================
#     Set up Problems
# ======================================================================
ap = AeroProblem(name=name+'_mach'+str(mach)+'_alpha'+str(alpha),
  mach=mach,
  reynolds=reynolds,
  reynoldsLength = reynoldsLength,
  T = T ,
  alpha=alpha,
  areaRef=areaRef,
  chordRef=chordRef,
  xRef=0.25,yRef=0.0,zRef=0.0,
  evalFuncs=['cl','cd','cmz']
)
tp = TransiProblem(name=name+'_mach'+str(mach)+'_alpha'+str(alpha),
  mach=mach,
  reynolds=reynolds/reynoldsLength,
  T=T,
  nCritTS=NTS,
  nCritCF=NCF,
  spanDirection=spanDirection,
  TurbulentIntencity = Tu*0.01,
  )
atp = AeroTransiProblem(ap,tp)
# ======================================================================
#     Set up Solvers
# ======================================================================
# if flags[aeroID]:
if 1==1:
  # we first try to not use flag
  print("OKay . Before CFDSOlver")
  CFDSolver = ADFLOW(options=aeroOptions,comm=comm)
  print("OKay . Inside CFDSOlver, and gcomm is",gcomm.rank,",comm is,",comm.rank)
  # CFDSolver.setDVGeo(DVGeo)
  # mesh = USMesh(options=meshOptions, comm=comm)
  # CFDSolver.setMesh(mesh)
  # pos = np.array([-0.4,-0.3,-0.25,-0.025,0.2,0.4])
  pos = np.array([0.5])
  CFDSolver.addSlices(spanDirection,pos,sliceType='absolute')
  transiSolver = None
if flags[transiID]:
  transiSolver = TransitionCalc(options=transiOptions,comm=comm)
  # CFDSolver = None
# define Aero transi solver
AT = AeroTransi(CFDSolver, transiSolver,gcomm,options=lfOptions)
# ======================================================================
#     Solve Functions:
# ======================================================================
funcs = {}
CFDSolver(ap)
CFDSolver.evalFunctions(ap,funcs)
# AT(atp)
# AT.evalFunctions(atp, funcs)
if(gcomm.rank == 0):
  print(funcs)