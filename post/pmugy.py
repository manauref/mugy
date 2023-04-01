#[ ................................................................ ]#
#[
#[ pmugy
#[ A collection of methods to be used in other scripts for
#[ post-processing pmugy data.
#[
#[ Some basic editing rules:
#[   - Use 2 spaces for indenting (no tab preferably).
#[
#[ ................................................................ ]#


from mpi4py import MPI
import numpy as np
import adios2
import sys

class pmIO:
  def __init__(self, **kwargs):
    self.ad = adios2.ADIOS()
    self.ioCounter = 0
    pass

  #[ Create an IO object:
  def declareIO(**kwargs):
    if 'name' in kwargs:
      return kwargs['name'], self.ad.DeclareIO(kwargs['name'])
    else:
      ob_name = 'ad_ioOb'+str(self.ioCounter) 
      self.ioCounter = self.ioCounter+1
      return ob_name, self.ad.DeclareIO(ob_name)
  #[ .......... end of declareIO method ........... ]#

  #[ Open a data file and return a handle to it.
  def fOpen(self, fileName, **kwargs):
    if 'ioObject' in kwargs:
      return _, kwargs['ioObject'],  kwargs['ioObject'].Open(fileName, adios2.Mode.Read)
    else:
      reader_name = kwargs.get('readerName', 'ad_reader'+str(self.ioCounter)) 
      self.ioCounter = self.ioCounter+1
      ad_read = self.ad.DeclareIO(reader_name)
      return reader_name, ad_read, ad_read.Open(fileName, adios2.Mode.Read)
  #[ .......... end of fOpen method ........... ]#

  #[ Close a file/stream.
  def fClose(self, stream):
    stream.Close()
  #[ .......... end of fClose method ........... ]#

  #[ Remove IO object:
  def removeIOobject(self, ioObName):
    return self.ad.RemoveIO(ioObName)
  #[ .......... end of removeIOobject method ........... ]#

  #[ Inquire a variable.
  def varInquire(self, varName):
    return self.ad_read.InquireVariable(varName)
  #[ .......... end of varInquire method ...... ]#

  #[ Get the shape of a variable.
  def varShape(self, varName, **kwargs):
    shapeOut = 0
    if 'fileName' in kwargs:
      ad_readNm, ad_read, ad_istream = self.fOpen(kwargs['fileName'])
      ad_varin = ad_read.InquireVariable(varName)
      shapeOut = ad_varin.Shape()
      self.fClose(ad_istream)
      closed = self.ad.RemoveIO(ad_readNm)
    elif 'ioObject' in kwargs:
      ad_varin = kwargs['ioObject'].InquireVariable(varName)
      shapeOut = ad_varin.Shape()
    else:
      sys.exit("varShape requires 'fileName' or 'ioObject' input.")
    return shapeOut
  #[ .......... end of varShape method ...... ]#

  #[ Translate an adios datatype to a numpy datatype.
  def adType2npType(self, adtype):
    if adtype == 'float':
      return np.single
    elif adtype == 'float complex':
      return np.csingle
    elif adtype == 'double':
      return np.double
    elif adtype == 'double complex':
      return np.cdouble
  #[ .......... end of adTyope2npType ....... ]#

  #[ Get the data type of a variable.
  def varType(self, varName, **kwargs):
    typeOut = ''
    if 'fileName' in kwargs:
      ad_readNm, ad_read, ad_istream = self.fOpen(kwargs['fileName'])
      ad_varin = ad_read.InquireVariable(varName)
      typeOut = ad_varin.Type()
      self.fClose(ad_istream)
      closed = self.ad.RemoveIO(ad_readNm)
    elif 'ioObject' in kwargs:
      ad_varin = kwargs['ioObject'].InquireVariable(varName)
      typeOut = ad_varin.Type()
    else:
      sys.exit("varType requires 'fileName' or 'ioObject' input.")

    numpyOut = kwargs.get('numpy', False)
    if numpyOut:
      return self.adType2npType(typeOut)
    else:
      return typeOut
  #[ .......... end of varShape method ...... ]#

  #[ Read a variable.
  def varRead(self, varName, **kwargs):
    #[ Establish file handle.
    closeFile = False
    if 'fileName' in kwargs:
      ad_readNm, ad_read, ad_iStream = self.fOpen(kwargs['fileName'])
      closeFile = True
    elif ('ioObject' in kwargs) and ('fileStream' in kwargs):
      ad_read    = kwargs['ioObject']
      ad_iStream = kwargs['fileStream']
    else:
      sys.exit("varRead requires 'fileName' or ('ioObject' and 'fileStream') input.")

    #[ Establish data selection shape.
    ad_varin = ad_read.InquireVariable(varName)
    if 'select' in kwargs:
      #[ 'select' argument must be a two-list list of starts and counts.
      ad_varin.SetSelection([kwargs['select'][0], kwargs['select'][1]])
    else:
      ad_varshape = ad_varin.Shape()
      ad_varin.SetSelection([[0 for i in range(len(ad_varshape))], ad_varshape])

    #[ Read variable.
    if 'array' in kwargs:
      ad_iStream.Get(ad_varin, kwargs['array'], adios2.Mode.Sync)
      if closeFile:
        self.fClose(ad_iStream)
        self.removeIOobject(ad_readNm)
    else:
      inSize = ad_varin.SelectionSize()
      var = np.zeros(inSize, dtype=self.adType2npType(ad_varin.Type()))
      ad_iStream.Get(ad_varin, var, adios2.Mode.Sync)
      if closeFile:
        self.fClose(ad_iStream)
        self.removeIOobject(ad_readNm)
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
