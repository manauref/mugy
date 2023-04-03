#[ ......................................................... #]
#[
#[ Plot a frame of a pmugy dataset.
#[
#[
#[ ......................................................... #]

import pmugy as pmc
import numpy as np
import matplotlib.pyplot as plt

#dataDir = '/Users/manaure/Documents/multiscale/code/mugy/src/'
dataDir = '/home/manaurer/multiscale/code/mugy/src/'
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

print(varShape, varType)
#[ Read variable in.
fldIn = np.zeros([varShape[2],varShape[3]], dtype=varType)
pm.varRead(varName, fileName=fileRoot, select=varSelect, array=fldIn)
print(fldIn)

#[ Assuming the data is in Fourier space, plot the square amplitude.
kx = [kxMin[0]*np.arange(varShape[2]), kxMin[1]*np.arange(varShape[3]), kxMin[2]*np.arange(varShape[1])]
kxNodal = [ np.append([kx[i][0]-kxMin[i]/2], 
            (0.5*(kx[i][:-1]+kx[i][1:])).tolist() + [kx[i][0]-kxMin[i]/2]) for i in range(len(kxMin))]
X = [np.outer(kxNodal[0],np.ones(np.size(kxNodal[1]))),
     np.outer(np.ones(np.size(kxNodal[0])),kxNodal[1])]

plt.pcolormesh(X[0], X[1], np.abs(fldIn))
plt.colorbar()
plt.show()

#[ Below we test opening and reading two files at the same time.

ad_readerNm, ad_reader, ad_stream = pm.fOpen(fileRoot)  #[ Open file.
ad_readerNm_phi, ad_reader_phi, ad_stream_phi = pm.fOpen(dataDir + 'phik.bp')  #[ Open file.

dataShape = pm.varShape(varName, ioObject=ad_reader)  #[ Variable shape.
#dataShape = pm.varShape(varName, fileName=fileRoot, readerName='read3')
#dataShape = pm.varShape(varName, fileName=fileRoot)
print(dataShape)
#
#dataShape = pm.varShape(varName, fileName=dataDir + 'phik.bp')  #[ Variable shape.
dataShape = pm.varShape(varName, ioObject=ad_reader_phi)  #[ Variable shape.
print(dataShape)

pm.fClose(ad_stream_phi)  #[ Close file.
pm.removeIOobject(ad_readerNm_phi)
pm.fClose(ad_stream)  #[ Close file.
pm.removeIOobject(ad_readerNm)
