#[ ................................................................ ]#
#[
#[ pmugy
#[ A collection of methods to be used in other scripts for
#[ post-processing pmugy simulations. 
#[
#[ Some basic editing rules:
#[     - Use 4 spaces for indenting (no tab preferably).
#[
#[ ................................................................ ]#


from mpi4py import MPI
import numpy as np
import adios2
import sys

class pmIO:
   def __init__(self, **kwargs):
       self.ad = adios2.ADIOS()
       self.ad_read = self.ad.DeclareIO(kwargs.get('readerName', 'ad_reader'))
       pass

   #[ Open a data file and return a handle to it.
   def fOpen(self, fileName):
       return self.ad_read.Open(fileName, adios2.Mode.Read)
   #[ .......... end of fOpen method ........... ]#

   #[ Close a file/stream.
   def fClose(self, ad_stream):
        ad_stream.Close()
   #[ .......... end of fOpen method ........... ]#

   #[ Inquire a variable.
   def varInquire(self, varName):
       return self.ad_read.InquireVariable(varName)
   #[ .......... end of varInquire method ...... ]#

   #[ Get the shape of a variable.
   def varShape(self, varName, **kwargs):
       if 'fileName' in kwargs:
           ad_istream = self.fOpen(kwargs['fileName'])
           ad_varin = self.ad_read.InquireVariable(varName)
           self.fClose(ad_istream)
       else:
           ad_varin = self.ad_read.InquireVariable(varName)
       return ad_varin.Shape()
   #[ .......... end of varShape method ...... ]#

   #[ Read a variable.
   def varRead(self, varName, **kwargs):
       #[ Establish file handle.
       closeFile = False
       if 'fileName' in kwargs:
           ad_iStream = self.ad_read.Open(kwargs['fileName'], adios2.Mode.Read)
           closeFile = True
       elif 'fileHandle' in kwargs:
           ad_iStream = kwargs['fileHandle']
       else:
           sys.exit("varRead requires 'fileName' or 'fileHandle' input.")

       #[ Establish data selection shape.
       ad_varin = self.ad_read.InquireVariable(varName)
       if 'select' in kwargs:
           #[ 'select' argument must be a two-list list of starts and counts.
           ad_varin.SetSelection([kwargs['select'][0], kwargs['select'][1]])
       else:
           ad_varshape = ad_varin.Shape()
           ad_varin.SetSelection([[0 for i in range(len(ad_varshape))], ad_varshape])

       #[ Read variable.
       if 'var' in kwargs:
           ad_iStream.Get(ad_varin, kwargs['var'], adios2.Mode.Sync)
           if closeFile:
               self.fClose(ad_iStream)
       else:
           inSize = ad_varin.SelectionSize()
           var = np.zeros(inSize, dtype=np.single)
           ad_iStream.Get(ad_varin, var, adios2.Mode.Sync)
           if closeFile:
               self.fClose(ad_iStream)
           return var
   #[ .......... end of varRead .............. ]#

#dataDir = '/Users/manaure/Documents/multiscale/code/mugy/src/'
#fileName = 'momk.bp'
#
##  debug mode
#adios = adios2.ADIOS(adios2.DebugON)
#ioRead = adios.DeclareIO("ioReader")
#
## ADIOS Engine
#ibpStream = ioRead.Open(dataDir + fileName, adios2.Mode.Read)
#
#advarin = ioRead.InquireVariable("momk")
#if advarin is not None:
#    advar_shape = advarin.Shape()
#    print(advar_shape)
#    advarin.SetSelection([[0 for i in range(4)], advar_shape])
#    inSize = advarin.SelectionSize()
#    print('Incoming size ' + str(inSize))
#
#    var_k = np.zeros(inSize, dtype=np.single)
#
#    ibpStream.Get(advarin, var_k, adios2.Mode.Sync)
#
##    print('Incoming temperature map')
##
##    for i in range(0, inTemperatures.size):
##        print(str(inTemperatures[i]) + ' ')
##
##        if (i + 1) % 4 == 0:
##            print()
#
#ibpStream.Close()
