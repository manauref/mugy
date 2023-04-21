#[ ......................................................... #]
#[ mugy postprocessing script:
#[
#[ Plot a frame of a fourier array, either in Fourier space or
#[ in real space by FFTing it, from a file that contains multiple
#[ frames.
#[
#[ Note that the example below may not be ideal because it
#[ opens and closes the file in every iteration of the loop.
#[ It'd be better to open the file once outside of the loop,
#[ do the reads inside the loop, and close it after the loop.
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

dataDir  = '/home/manaurer/multiscale/code/mugy/src/'
fileName = 'phik.bp'
varName  = 'globalVariable'

space = 'real'  #[ Choose 'real' or 'fourier' space.

#[ Indicate the first frame (frameSelect[0])
#[ and the number of frames (frameSelect[1]) desired.
frameSelect = [0, 1]


#[ .......... End of user inputs (MAYBE) ............. ]#

filePathName = dataDir + fileName

pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ Get the grid (note that pcolormesh expects it to be nodal).
if space == 'real':
  xNodal = pm.xGrid(fileName=filePathName, nodal=True)
else:
  kxNodal = pm.kGrid(fileName=filePathName, nodal=True)
  xNodal, kx0centerRollNum = pm.kx0centerGrid(kxNodal)
  
X = [np.outer(xNodal[0],np.ones(np.size(xNodal[1]))),
     np.outer(np.ones(np.size(xNodal[0])),xNodal[1])]

#[ Shape and datatype of the variable.
varShape = pm.varShape(varName, fileName=filePathName)
varType  = pm.varType(varName, fileName=filePathName, numpy=True)

#[ Assume the variable have the shape [Nkz x Nkx x Nky].
#[ Select one x-y plane with dimSelect = [starts, counts].
dimSelect = [[0,0,0], [1,varShape[1],varShape[2]]]

#[ Preallocate space for a single frame of the variable.
fldIn = np.zeros([1, varShape[1],varShape[2]], dtype=varType)

for frIdx in range(frameSelect[0], frameSelect[0]+frameSelect[1]):

  #[ Read variable in.
  stepSelect = [frIdx, 1]
  pm.varReadSteps(varName, fileName=filePathName, select=dimSelect, selectSteps=stepSelect, array=fldIn)
  print("frIdx = ",frIdx)
  print(fldIn[0,0,:])
  print("")

  if space == 'real':
    #[ Transform variable to real space.
    fldOut = np.fft.irfft2(np.squeeze(fldIn), s=[np.size(xNodal[0])-1, np.size(xNodal[1])-1])
  else:
    #[ Keep in fourier space, but reorganize so kx's are in ascending order.
    fldOut = pm.kx0centerVar(np.squeeze(fldIn), rollnum=kx0centerRollNum)
    fldOut = np.abs(fldOut)  #[ Plot the complex amplitude.
  
  #[ Create plot
  plt.pcolormesh(X[0], X[1], fldOut)
  plt.colorbar()
  plt.show()

