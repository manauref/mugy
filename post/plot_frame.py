#[ ......................................................... #]
#[
#[ Plot a frame of a pmugy dataset.
#[
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

dataDir = '/Users/manaure/Documents/multiscale/code/mugy/src/'
fileName = 'momk.bp'
#fileName = 'phik.bp'

varName = 'momk'

kxMin = [0.12566, 0.12566, 1.0]


#[ .......... End of user inputs (MAYBE) ............. ]#

fileRoot = dataDir + fileName


pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ Shape and datatype of the variable.
varShape = pm.varShape(varName, fileName=fileRoot)
varType  = pm.varType(varName, fileName=fileRoot, numpy=True)

#[ Moment variables have the shape [num moments x Nkz x Nkx x Nky ].
#[ Select one x-y plane:
varSelect = [[0,0,0,0], [1,1,varShape[2],varShape[3]]]

#[ Read variable in.
var = np.zeros([varShape[2],varShape[3]], dtype=varType)
pm.varRead(varName, fileName=fileRoot, select=varSelect, array=var)

#[ Assuming the data is in Fourier space, plot the square amplitude.
kx = [kxMin[0]*np.arange(varShape[2]), kxMin[1]*np.arange(varShape[3]), kxMin[2]*np.arange(varShape[1])]
kxNodal = [ np.append([kx[i][0]-kxMin[i]/2], 
            (0.5*(kx[i][:-1]+kx[i][1:])).tolist() + [kx[i][0]-kxMin[i]/2]) for i in range(len(kxMin))]
X = [np.outer(kxNodal[0],np.ones(np.size(kxNodal[1]))),
     np.outer(np.ones(np.size(kxNodal[0])),kxNodal[1])]

plt.pcolormesh(X[0], X[1], np.abs(var))
plt.colorbar()
plt.show()

#ad_stream = pm.fOpen(fileRoot)  #[ Open file.
#ad_stream_phi = pm.fOpen(dataDir + 'phik.bp')  #[ Open file.

#dataShape = pm.varShape(varName)  #[ Variable shape.

#dataShape = pm.varShape(varName, fileName=dataDir + 'phik.bp')  #[ Variable shape.
#print(dataShape)

#pm.fClose(ad_stream_phi)  #[ Close file.
#pm.fClose(ad_stream)  #[ Close file.
