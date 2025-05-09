# This is a template that should be used for setting up
# RANS analysis scripts

# ======================================================================
#         Import modules
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

parser = argparse.ArgumentParser()
parser.add_argument("--Chord", type=float, default=1.1)
parser.add_argument("--Le", type=float, default=0.5)
parser.add_argument("--Te", type=float, default=0.5)
parser.add_argument("--mach", type=str, default=0.19)
# parser.add_argument("--mach", type=str, default=2.0)
parser.add_argument("--procs", help="number of processors", type=int, default=16)
parser.add_argument('--alpha', type=str, default=0.7498)
parser.add_argument("--reynolds", type=float, default=5.6e6)
args = parser.parse_args()

# ======================================================================
#         Input Information -- Modify accordingly!
# ======================================================================
# set the aerodynamic data
outputDirectory = './output/'
try:
    os.stat(outputDirectory)
except:
    os.mkdir(outputDirectory) 
iter = 0
# set the aerodynamic data
alpha = args.alpha
mach = args.mach
reynolds = args.reynolds
name = 'test'
T = 288
 
spanDirection = 'z'
reynoldsLength = args.Chord

chord = args.Chord
alpha_Le = args.Le
alpha_Te = args.Te


FileName_Input = ['slice.dat']
FileName = FileName_Input[0]

#Get the slice data from transi_run_data_slices.txt
Xdata1,Ydata1,Zdata1,CPdata1 = numpy.loadtxt(FileName,unpack=True, usecols = (0,1,2,12))


# we are hardcoding these here, but in the actual code these will be provided by the shape
# of the array coming out of ADflow or whatever other source provides the data
aeroData = {}
aeroData['Slice_0001'] = {}
aeroData['Slice_0001']['data'] = []
aeroData['Slice_0001']['data'].append({})
aeroData['Slice_0001']['data'][0]['x'] = Xdata1
aeroData['Slice_0001']['data'][0]['y'] = Ydata1
aeroData['Slice_0001']['data'][0]['z'] = Zdata1
aeroData['Slice_0001']['data'][0]['CoefPressure'] = CPdata1.copy()
aeroData['Slice_0001']['data'][0]['chord'] = []
aeroData['Slice_0001']['data'][0]['chord'] = chord

# FileName = FileName_Input[1]
# Xdata2,Ydata2,Zdata2,CPdata2 = numpy.loadtxt(FileName,unpack=True, usecols = (0,1,2,12))

# aeroData['Slice_0002'] = {}
# aeroData['Slice_0002']['data'] = []
# aeroData['Slice_0002']['data'].append({})
# aeroData['Slice_0002']['data'][0]['x'] = Xdata2
# aeroData['Slice_0002']['data'][0]['y'] = Ydata2
# aeroData['Slice_0002']['data'][0]['z'] = Zdata2
# aeroData['Slice_0002']['data'][0]['CoefPressure'] = CPdata1.copy()
# aeroData['Slice_0002']['data'][0]['chord'] = []
# aeroData['Slice_0002']['data'][0]['chord'] = 1.0
# trLoaOld = numpy.array([[[0.25779236,  0.05959119,  0.5       ,  0.05      ],
#                          [0.39905682, -0.05806271,  0.5       ,  0.05  ]    ,
#                          [0.39905682, -0.05806271,  0.5       ,  0.05  ]    ,
#                          [0.39905682, -0.05806271,  0.5       ,  0.05  ]     ]])
# trLoaIter = trLoaOld
# nCritTS - Critical Nfactor for TS Wave
# nCritCF - Critical Nfactor for CF Wave
# Both computed from surface roughness and turbulent intensity using experimental data
# We need to compute phiLE and phiTE from the geometry in ADFlow so that they will be
# affected properly by the geometric design variables..
# these should be a varable of size [nPart,nSpan]# Also if we actually need the reynolds number, this should inherit from aeroproblem
tr = TransiProblem(name=name, mach=mach, reynolds=reynolds/reynoldsLength, T=T,
                  nCritTS=6,nCritCF=5,TurbulentIntencity = 0.0007, spanDirection=spanDirection)

transiOptions = {
  #Relaxation Parameter
  'outputDirectory':outputDirectory,
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

Transi = TransitionCalc(options = transiOptions)
Transi.AeroToTransi(aeroData)

# Main call does everything
Transi(tr)

TrLoA, TrDir, xMinMax = Transi.getADflowValue()
print('Transition Data:',TrLoA,TrDir,xMinMax)

nSlice = 1
for i in range(1):
    for j in range(2*nSlice):
        if(j<nSlice):
            k = j
        else:
            k = nSlice-1 - j
        x = (TrLoA[i][j][0] - xMinMax[i][k][1])/(xMinMax[i][k][0]-xMinMax[i][k][1])
        print("Normalized transition locations",x)
print(TrLoA[0][0][0],TrLoA[0][1][0])

