#[ ......................................................... #]
#[
#[ Plot a frame of a pmugy dataset.
#[
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

dataDir  = '/home/manaurer/multiscale/code/mugy/src/'
#fileName = 'momk.bp'
#varName  = 'momk'
fileName = 'arr.bp'
varName  = 'globalVariable'


#[ .......... End of user inputs (MAYBE) ............. ]#

filePathName = dataDir + fileName

pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ Shape and datatype of the variable.
varShape = pm.varShape(varName, fileName=filePathName)
varType  = pm.varType(varName, fileName=filePathName, numpy=True)

#[ Moment variables have the shape [numMoments x Nkz x Nkx x Nky ],
#[ where numMoments is the total number of moments across all species).
#[ Select one x-y plane with varSelect = [starts, counts].
#varSelect = [[0,0,0,0], [1,1,varShape[2],varShape[3]]]
#
##[ Read variable in.
#fldIn = np.zeros([varShape[2],varShape[3]], dtype=varType)
#pm.varRead(varName, fileName=filePathName, select=varSelect, array=fldIn)
#
##[ Assuming the data is in Fourier space, plot the square amplitude.
#kxNodal = pm.kGrid(fileName=filePathName, nodal=True)
#X = [np.outer(kxNodal[0],np.ones(np.size(kxNodal[1]))),
#     np.outer(np.ones(np.size(kxNodal[0])),kxNodal[1])]
#
##[ Create plot
#plt.pcolormesh(X[0], X[1], np.abs(fldIn))
#plt.colorbar()
#plt.show()


##[ Real array.
varSelect = [[0,0,0], [1,varShape[1],varShape[2]]]
fldIn = np.zeros([varShape[1],varShape[2]], dtype=varType)
pm.varRead(varName, fileName=filePathName, select=varSelect, array=fldIn)
xNodal = pm.xGrid(fileName=filePathName, nodal=True)
X = [np.outer(xNodal[0],np.ones(np.size(xNodal[1]))),
     np.outer(np.ones(np.size(xNodal[0])),xNodal[1])]
print(fldIn[:,21])
plt.pcolormesh(X[0], X[1], fldIn)
plt.colorbar()
plt.show()

# Test fft on 43x22 grid
gld = np.zeros([varShape[1],varShape[2]], dtype=varType)
kxMin = pm.attrRead('kxMin',fileName=dataDir+'momk.bp')
xC = pm.xGrid(fileName=filePathName)
for i in range(varShape[1]):
  for j in range(varShape[2]):
    gld[i,j] += 2.5e-2*np.sin(kxMin[0]*xC[0][i])*np.cos(kxMin[1]*xC[1][j]);

#print("gld[:,21] = ",gld[:,21])
#gldk = np.fft.rfft2(gld)
#print(gldk)
#gld = np.fft.irfft2(gldk)
##print("gld[:,21] = ",gld[:,21])

#plt.pcolormesh(X[0], X[1], fldIn-gld)
#plt.colorbar()
#plt.show()

#print("Lx = ",xC[0][-1]-xC[0][0],xC[1][-1]-xC[1][0],xC[2][-1]-xC[2][0])
#print("x[0] = ",xC[0])
#print("x[1] = ",xC[1])
