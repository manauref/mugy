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


#[ .......... End of user inputs (MAYBE) ............. ]#

fileRoot = dataDir + fileName


pm = pmc.pmIO()  #[ Initialize the pmugy class

#[ Shape of the variable.
varShape = pm.varShape(varName, fileName=fileRoot)  #[ Variable shape.

#[ Moment variables have the shape [num moments x Nkz x Nkx x Nky ].
#[ Select one x-y plane:
varSelect = [[0,0,0,0], [1,1,varShape[2],varShape[3]]]

#[ Read variable in.
var = np.zeros([varShape[2],varShape[3]])
pm.varRead(varName, fileName=fileRoot, select=varSelect, var=var)

X = [np.outer(np.arange(varShape[2]),np.ones(varShape[3])),
     np.outer(np.ones(varShape[2]),np.arange(varShape[3]))]

print(var)
plt.pcolormesh(X[0], X[1], var)
plt.show()

#ad_stream = pm.fOpen(fileRoot)  #[ Open file.
#ad_stream_phi = pm.fOpen(dataDir + 'phik.bp')  #[ Open file.

#dataShape = pm.varShape(varName)  #[ Variable shape.

#dataShape = pm.varShape(varName, fileName=dataDir + 'phik.bp')  #[ Variable shape.
#print(dataShape)

#pm.fClose(ad_stream_phi)  #[ Close file.
#pm.fClose(ad_stream)  #[ Close file.
