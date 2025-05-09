import os
import numpy as np
from pathlib import Path
from collections import OrderedDict
import re


# Tur slices
# print("Where am I", flush=True)
os.makedirs('./Slice_pure',exist_ok=True)

fi = open('./output_aero/test_mach0.19_alpha0.7498_000_slices.dat','r')
lines = fi.readlines() # have read all the lines

Slice = [0.50000000]
for m in range(len(Slice)):
    n = Slice[m]
    for j,line in enumerate(lines):
        if '"Slice_0001 None Absolute Regular z = '+str(n) in line:
        # if 'Zone T= ' in line:
            flagline = j
            # Parts = lines[j+1].split(',')
            # nNodes = int(Parts[0][7:].replace(',',' '))
            ntail  = int(lines[j+1].split()[4])
            nNodes = int(lines[j+1].split()[2])
        
    print('-> nNodes:',nNodes,'ntail:',ntail)
    # print('-> nNodes:',nNodes, flush=True)
    nVars=18

    SliceSurfData= np.zeros((nNodes,nVars))
    # tail= np.zeros((ntail,2))
    num = 0
    for k in range(flagline+3,flagline+nNodes+3):
        if 'E-' in lines[k] or 'E+' in lines[k]:
            SliceSurfData[num,:] = np.array(list(map(float, lines[k].split())))
            num += 1
            # tailbeg = k+1
    # tmp = 0
    # for k in range(flagline+nNodes+3,flagline+nNodes+3+ntail):
    #     tail[tmp,:] = np.array(list(map(int, lines[k].split())))
    #     tmp += 1
    # print('-> dataline:',num,' tailline:',tmp)
    np.savetxt('./Slice_pure/slice_'+str(m+1)+'.dat', SliceSurfData)

# Transi
# casename = 'Transi'
# os.chdir(casename)
chordlength = [1.0]
alpha_le = [0.0]
alpha_te = [0.0]

for i in range(len(chordlength)): # len(chordlength)
    # print(f"now i = {str(i)}", flush=True)
    Fre = [0]
    casename = 'Slice_'+str(i+1)
    os.makedirs(casename,exist_ok=True)
    os.chdir(casename)
    for m in range(len(Fre)):
        # print(f"now m = {str(m)}", flush=True)
        casename = 'Fre'+str(Fre[m])
        os.makedirs(casename,exist_ok=True)
        os.chdir(casename)
        os.system('cp ../../3d.py 3d.py')
        os.system('cp ../../Slice_pure/slice_'+str(i+1)+'.dat'+' slice.dat')
        os.system('python 3d.py  --Chord {0}'.format(chordlength[i])+' --Le {0}'.format(alpha_le[i])+' --Te {0}'.format(alpha_te[i])+f' 2>&1|tee print{str(i)}_{str(m)}.log')
        os.chdir(os.path.pardir)
    os.chdir(os.path.pardir)

#transiLoc
Error = [100]
Fre = [0]
os.makedirs('TransiLoc',exist_ok=True)
SliceNum = 1
Loc_u = np.zeros((SliceNum,4))
Loc_l = np.zeros((SliceNum,4))
for i in range(SliceNum):
    transiLoc_upper = np.zeros([len(Fre),4])
    transiLoc_lower = np.zeros([len(Fre),4])
    casename = 'Slice_'+str(i+1)
    os.chdir(casename)
    for m in range(len(Fre)):
        casename = 'Fre'+str(Fre[m])
        os.chdir(casename)
        os.chdir('output')
        myfile = Path('./transiLoc.dat')
        if myfile.exists():
            transiLoc = np.loadtxt('./transiLoc.dat')
            transiLoc_upper[m,:] = transiLoc[0,:]
            transiLoc_lower[m,:] = transiLoc[1,:]
            os.chdir(os.path.pardir)
            os.chdir(os.path.pardir)
        else:
            transiLoc_upper[m,:] = Error[:]
            transiLoc_lower[m,:] = Error[:]
            os.chdir(os.path.pardir)
            os.chdir(os.path.pardir)
    # os.chdir(os.path.pardir)
    np.savetxt('../TransiLoc/transiLoc_upper_slice'+str(i+1)+'.dat',transiLoc_upper)
    np.savetxt('../TransiLoc/transiLoc_lower_slice'+str(i+1)+'.dat',transiLoc_lower)

    index_upper=np.argmin(transiLoc_upper, axis=0)[0]
    index_lower=np.argmin(transiLoc_lower, axis=0)[0]

    Loc_upper = transiLoc_upper[index_upper,:]
    Loc_lower = transiLoc_lower[index_lower,:]
    np.savetxt('../TransiLoc/Loc_upper_slice_'+str(i+1)+'.dat',Loc_upper)
    np.savetxt('../TransiLoc/Loc_lower_slice_'+str(i+1)+'.dat',Loc_lower)

    Loc_u[i,:] = Loc_upper
    Loc_l[SliceNum-i-1,:] = Loc_lower
    os.chdir(os.path.pardir)

Loc = np.vstack((Loc_u,Loc_l))
np.savetxt('./TransiLoc/Loc.dat',Loc)







