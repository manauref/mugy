#[ ......................................................... #]
#[
#[ Plot a fourier array, either in Fourier space or in
#[ real space by FFTing it.
#[
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

dataDir  = '/home/manaurer/multiscale/code/mugy/src/'
fileName = 'phik.bp'
varName  = 'globalVariable'

space = 'real'  #[ Choose 'real' or 'fourier' space.


#[ .......... End of user inputs (MAYBE) ............. ]#

filePathName = dataDir + fileName

pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ Shape and datatype of the variable.
varShape = pm.varShape(varName, fileName=filePathName)
varType  = pm.varType(varName, fileName=filePathName, numpy=True)

#[ Assume the variable have the shape [Nkz x Nkx x Nky].
#[ Select one x-y plane with varSelect = [starts, counts].
varSelect = [[0,0,0], [1,varShape[1],varShape[2]]]

#[ Read variable in.
fldIn = np.zeros([varShape[1],varShape[2]], dtype=varType)
pm.varRead(varName, fileName=filePathName, select=varSelect, array=fldIn)

#[ Assuming the data is in Fourier space, plot the square amplitude.
if space == 'real':
  #[ Transform variable to real space and 
  xNodal = pm.xGrid(fileName=filePathName, nodal=True)
  fldOut = np.fft.irfft2(fldIn, s=[np.size(xNodal[0])-1, np.size(xNodal[1])-1])
else:
  #[ Keep in fourier space, but reorganize so kx's are in ascending order.
  kxNodal = pm.kGrid(fileName=filePathName, nodal=True)
  xNodal, fldOut = pm.kx0center(kxNodal, fldIn)
  fldOut = np.abs(fldOut)

X = [np.outer(xNodal[0],np.ones(np.size(xNodal[1]))),
     np.outer(np.ones(np.size(xNodal[0])),xNodal[1])]

#[ Create plot
plt.pcolormesh(X[0], X[1], fldOut)
plt.colorbar()
plt.show()

