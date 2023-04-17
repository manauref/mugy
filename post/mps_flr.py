#[ ......................................................... #]
#[ mugy postprocessing script
#[
#[ Plot the FLR functions computed with mugy in C and/or Python.
#[
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

dataDir  = '/home/manaurer/multiscale/code/mugy/src/'
fileName = ['pbFLRop.bp', 'poissonDb.bp', 'poissonSb.bp']
varName  = 'globalVariable'

plotSpec = "ion"

muMass = {"elc" : 60.59, "ion" : 1.0}
tau    = {"elc" : 1.0  , "ion" : 1.0}

#[ .......... End of user inputs (MAYBE) ............. ]#

xyLabelFontSize = 17
legendFontSize  = 14

pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ ............ Poisson bracket FLR operators ........... ]#
filePathName = [dataDir+fileName[0],dataDir+fileName[1],dataDir+fileName[2],]

#[ Shape and datatype of the variable.
varShape = pm.varShape(varName, fileName=filePathName[0])
varType  = pm.varType(varName, fileName=filePathName[0], numpy=True)

#[ The quantity pbFLRop has the FLR operators that act on variables
#[ going into Poisson brackets. These are
#[   <J_0>=Gamma_0^{1/2}, 0.5*hatLap <J_0>, (1+0.5*hatLap+hathatLap) <J_0>.
#[ For each species, so the array has shape [(numSpecies x 3) x Nkx x Nky ].

specOff = 0
if plotSpec == "ion":
  specOff = 3

#[ Read <J_0>.
varSelect = [[specOff+0,0,0], [1,varShape[1],varShape[2]]]
pbFLR0c = np.zeros([varShape[1],varShape[2]], dtype=varType)
pm.varRead(varName, fileName=filePathName[0], select=varSelect, array=pbFLR0c)
#[ Read 0.5*hatLap <J_0>.
varSelect = [[specOff+1,0,0], [1,varShape[1],varShape[2]]]
pbFLR1c = np.zeros([varShape[1],varShape[2]], dtype=varType)
pm.varRead(varName, fileName=filePathName[0], select=varSelect, array=pbFLR1c)
#[ Read (1+0.5*hatLap+hathatLap) <J_0>.
varSelect = [[specOff+2,0,0], [1,varShape[1],varShape[2]]]
pbFLR2c = np.zeros([varShape[1],varShape[2]], dtype=varType)
pm.varRead(varName, fileName=filePathName[0], select=varSelect, array=pbFLR2c)

#[ Plot FLR function at kx=0.
kx = pm.kGrid(fileName=filePathName[0], nodal=False)

#[ Compute FLR functions with python.
flrOps = pm.calcFLR(kx,tau[plotSpec],muMass[plotSpec],only=["avgJ0","hatLap","hathatLap"])

pbFLR0py = flrOps["avgJ0"]
pbFLR1py = 0.5*flrOps["hatLap"]*flrOps["avgJ0"]
pbFLR2py = (1.+0.5*flrOps["hatLap"]+flrOps["hathatLap"])*flrOps["avgJ0"]

#[ Prepare figure.
figProp = (6.4,6.8)
axPos   = [[0.18, 0.690, 0.725, 0.285],
           [0.18, 0.395, 0.725, 0.285],
           [0.18, 0.100, 0.725, 0.285]]
fig     = plt.figure(figsize=figProp)
ax      = [fig.add_axes(pos) for pos in axPos]

ax[0].plot(kx[1], pbFLR0c[0,:], linestyle='-')
ax[0].plot(kx[1], pbFLR0py[0,:], linestyle='--')
ax[1].plot(kx[1], pbFLR1c[0,:], linestyle='-')
ax[1].plot(kx[1], pbFLR1py[0,:], linestyle='--')
ax[2].plot(kx[1], pbFLR2c[0,:], linestyle='-')
ax[2].plot(kx[1], pbFLR2py[0,:], linestyle='--')

ax[2].set_xlabel(r'$k_y\rho_s$', fontsize=xyLabelFontSize)
ax[0].set_ylabel(r'$\left\langle J_0\right\rangle$', fontsize=xyLabelFontSize)
ax[1].set_ylabel(r'$\frac{1}{2}\hat{\nabla}_\perp\left\langle J_0\right\rangle$', fontsize=xyLabelFontSize)
ax[2].set_ylabel(r'$\left(1+\frac{1}{2}\hat{\nabla}_\perp+\hat{\hat{\nabla}}_\perp\right)\left\langle J_0\right\rangle$', fontsize=xyLabelFontSize)
ax[0].legend([r'mugy',r'Python'], fontsize=legendFontSize, frameon=False)
for i in range(2):
  plt.setp( ax[i].get_xticklabels(), visible=False)
plt.show()

##[ ............ Poisson equation FLR operators ........... ]#
#
#specOff = 0
#if plotSpec == "ion":
#  specOff = 1
#
##[ Read Db
#varShape = pm.varShape(varName, fileName=filePathName[1])
#varType  = pm.varType(varName, fileName=filePathName[1], numpy=True)
#varSelect = [[specOff+0,0,0], [1,varShape[1],varShape[2]]]
#Dbc = np.zeros([varShape[1],varShape[2]], dtype=varType)
#pm.varRead(varName, fileName=filePathName[1], select=varSelect, array=Dbc)
#
##[ Read Sb
#varShape = pm.varShape(varName, fileName=filePathName[2])
#varType  = pm.varType(varName, fileName=filePathName[2], numpy=True)
#varSelect = [[specOff+0,0,0], [1,varShape[1],varShape[2]]]
#Sbc = np.zeros([varShape[1],varShape[2]], dtype=varType)
#pm.varRead(varName, fileName=filePathName[2], select=varSelect, array=Sbc)
#
##[ Compute FLR functions with python.
#flrOps = pm.calcFLR(kx,tau[plotSpec],muMass[plotSpec],only=["Db","Sb"])
#
#Dbpy = flrOps["Db"]
#Sbpy = flrOps["Sb"]
#
##[ Prepare figure.
#figProp = (6.4,4.8)
#axPos   = [[0.18, 0.540, 0.725, 0.385],
#           [0.18, 0.140, 0.725, 0.385]]
#fig     = plt.figure(figsize=figProp)
#ax      = [fig.add_axes(pos) for pos in axPos]
#
#ax[0].plot(kx[1], Dbc[0,:], linestyle='-')
#ax[0].plot(kx[1], Dbpy[0,:], linestyle='--')
#ax[1].plot(kx[1], Sbc[0,:], linestyle='-')
#ax[1].plot(kx[1], Sbpy[0,:], linestyle='--')
#
#ax[1].set_xlabel(r'$k_y\rho_s$', fontsize=xyLabelFontSize)
#ax[0].set_ylabel(r'$S_b$', fontsize=xyLabelFontSize)
#ax[1].set_ylabel(r'$D_b$', fontsize=xyLabelFontSize)
#ax[0].legend([r'mugy',r'Python'], fontsize=legendFontSize, frameon=False)
#plt.setp( ax[0].get_xticklabels(), visible=False)
#plt.show()
