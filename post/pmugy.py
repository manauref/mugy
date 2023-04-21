#[ ................................................................ ]#
#[
#[ pmugy
#[ A collection of methods to be used in other scripts for
#[ post-processing pmugy data.
#[
#[ Some basic editing rules:
#[   - Use 2 spaces for indenting (no tab preferably).
#[
#[ Functions in this file
#[   declareIO
#[   fOpen
#[   fClose
#[   removeIOobject
#[   attrInquire
#[   attrRead
#[   attrHas
#[   varInquire
#[   varShape
#[   adType2npType
#[   varType
#[   varRead
#[   varReadSteps
#[   kGrid
#[   xGrid
#[   calcFLR
#[   kx0centerGrid
#[   kx0centerVar
#[
#[ ................................................................ ]#


from mpi4py import MPI
import numpy as np
import adios2
import sys
import scipy.special as scsp  #[ For Bessel functions.

class pmIO:
  def __init__(self, **kwargs):
    self.ad = adios2.ADIOS()
    self.ioCounter = 0 #[ Counts the IO engines created.
    pass

  #[ Create an IO object:
  def declareIO(**kwargs):
    if 'name' in kwargs:
      return kwargs['name'], self.ad.DeclareIO(kwargs['name'])
    else:
      ob_name = 'ad_io'+str(self.ioCounter) 
      self.ioCounter = self.ioCounter+1
      return ob_name, self.ad.DeclareIO(ob_name)
  #[ .......... end of declareIO method ........... ]#

  #[ Open a data file and return a handle to it.
  def fOpen(self, fileName, **kwargs):
    if 'io' in kwargs:
      #[ The return argument here is the name, but not sure if can extract from io.
      return _, kwargs['io'],  kwargs['io'].Open(fileName, adios2.Mode.Read)
    else:
      ad_io_name = kwargs.get('readerName', 'ad_reader_io'+str(self.ioCounter)) 
      self.ioCounter = self.ioCounter+1
      ad_io = self.ad.DeclareIO(ad_io_name)
      ad_eng = ad_io.Open(fileName, adios2.Mode.Read)
      return ad_io_name, ad_io, ad_eng
  #[ .......... end of fOpen method ........... ]#

  #[ Close a file/engine.
  def fClose(self, engine):
    engine.Close()
  #[ .......... end of fClose method ........... ]#

  #[ Remove IO object:
  def removeIOobject(self, ioObName):
    return self.ad.RemoveIO(ioObName)
  #[ .......... end of removeIOobject method ........... ]#

  #[ Get a dictionary of available attributes.
  def attrAvailable(self, **kwargs):
    ad_attr, openedFile = 0, False
    if 'fileName' in kwargs:
      ad_ioNm, ad_io, ad_eng = self.fOpen(kwargs['fileName'])
      ad_attrs = ad_io.AvailableAttributes()
    elif 'io' in kwargs:
      ad_attrs = kwargs['io'].AvailableAttributes()
    else:
      sys.exit("attrAvailable requires 'fileName' or 'io' input.")
    return ad_attrs, openedFile, ad_ioNm, ad_io, ad_eng
  #[ .......... end of attrAvailable .......... ]#

  #[ Inquire an attribute.
  def attrInquire(self, attrName, **kwargs):
    ad_attr, openedFile = 0, False
    if 'fileName' in kwargs:
      ad_ioNm, ad_io, ad_eng = self.fOpen(kwargs['fileName'])
      ad_attr = ad_io.InquireAttribute(attrName)
    elif 'io' in kwargs:
      ad_attr = kwargs['io'].InquireAttribute(attrName)
    else:
      sys.exit("attrInquire requires 'fileName' or 'io' input.")
    return ad_attr, openedFile, ad_ioNm, ad_io, ad_eng
  #[ .......... end of attrInquire method ...... ]#

  #[ Check if an attribute exists in file.
  def attrHas(self, attrName, **kwargs):
    ad_attrs, openedFile, ad_ioNm, ad_io, ad_eng = self.attrAvailable(**kwargs)
    hasAttr = attrName in ad_attrs.keys()
    if openedFile:
      self.fClose(ad_eng)
      closed = self.ad.RemoveIO(ad_ioNm)
    return hasAttr
  #[ .......... end of attrHas .......... ]#  

  #[ Read an attribute.
  def attrRead(self, attrName, **kwargs):
    #[ Establish file handle.
    ad_attr, openedFile, ad_ioNm, ad_io, ad_eng = self.attrInquire(attrName, **kwargs)
    if ad_attr.Type() == "string":
      val = ad_attr.DataString()[0]
    else:
      val = ad_attr.Data()
    if openedFile:
      self.fClose(ad_eng)
      closed = self.ad.RemoveIO(ad_ioNm)
    return val
  #[ .......... end of attrRead .............. ]#

  #[ Inquire a variable.
  def varInquire(self, varName, **kwargs):
    ad_var, openedFile = 0, False
    if 'fileName' in kwargs:
      ad_ioNm, ad_io, ad_eng = self.fOpen(kwargs['fileName'])
      ad_var = ad_io.InquireVariable(varName)
      openedFile = True
    elif 'io' in kwargs:
      ad_var = kwargs['io'].InquireVariable(varName)
    else:
      sys.exit("varInquire requires 'fileName' or 'io' input.")
    return ad_var, openedFile, ad_ioNm, ad_io, ad_eng
  #[ .......... end of varInquire method ...... ]#

  #[ Get the shape of a variable.
  def varShape(self, varName, **kwargs):
    ad_var, openedFile, ad_ioNm, ad_io, ad_eng = self.varInquire(varName, **kwargs)
    shape = ad_var.Shape()
    if openedFile:
      self.fClose(ad_eng)
      closed = self.ad.RemoveIO(ad_ioNm)
    return shape 
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
  #[ .......... end of adType2npType ....... ]#

  #[ Get the data type of a variable.
  def varType(self, varName, **kwargs):
    ad_var, openedFile, ad_ioNm, ad_io, ad_eng = self.varInquire(varName, **kwargs)
    typeOut  = ad_var.Type()
    if openedFile:
      self.fClose(ad_eng)
      closed = self.ad.RemoveIO(ad_ioNm)

    numpyOut = kwargs.get('numpy', False)
    if numpyOut:
      return self.adType2npType(typeOut)
    else:
      return typeOut
  #[ .......... end of varType method ...... ]#

  #[ Read a variable.
  def varRead(self, varName, **kwargs):
    if 'io' in kwargs:
      if not ('eng' in kwargs):
        sys.exit("varRead: is providing 'io' must also provide 'eng'.")
    ad_var, openedFile, ad_ioNm, ad_io, ad_eng = self.varInquire(varName, **kwargs)

    #[ Establish data selection shape.
    if 'select' in kwargs:
      #[ 'select' argument must be a two-list list of starts and counts.
      ad_var.SetSelection([kwargs['select'][0], kwargs['select'][1]])
      shape = kwargs['select'][1]
    else:
      shape = ad_var.Shape()
      ad_var.SetSelection([[0 for i in range(len(shape))], shape])

    #[ Read variable.
    var = kwargs.get('array',0)
    if not ('array' in kwargs):
      var = np.zeros(shape, dtype=self.adType2npType(ad_var.Type()))

    ad_eng.Get(ad_var, var, adios2.Mode.Sync)
    if openedFile:
      self.fClose(ad_eng)
      self.removeIOobject(ad_ioNm)
    return var
  #[ .......... end of varRead .............. ]#

  #[ Read a number of steps of a named variable.
  def varReadSteps(self, varName, **kwargs):
    if 'io' in kwargs:
      if not ('eng' in kwargs):
        sys.exit("varRead: is providing 'io' must also provide 'eng'.")
    ad_var, openedFile, ad_ioNm, ad_io, ad_eng = self.varInquire(varName, **kwargs)

    #[ Establish data selection shape.
    if 'select' in kwargs:
      #[ 'select' argument must be a two-list list of starts and counts.
      ad_var.SetSelection([kwargs['select'][0], kwargs['select'][1]])
      shape = kwargs['select'][1]
    else:
      shape = ad_var.Shape()
      ad_var.SetSelection([[0 for i in range(len(shape))], shape])
    #[ Select the desired steps.
    if 'selectSteps' in kwargs:
      #[ 'selectSteps' argument must be a two-list list of start and count.
      ad_var.SetStepSelection([kwargs['selectSteps'][0], kwargs['selectSteps'][1]])
    else:
      steps = ad_var.Steps()
      ad_var.SetStepSelection([0, steps])

    #[ Read variable.
    var = kwargs.get('array',0)
    if not ('array' in kwargs):
      var = np.zeros([steps]+shape, dtype=self.adType2npType(ad_var.Type()))

    ad_eng.Get(ad_var, var, adios2.Mode.Sync)
    if openedFile:
      self.fClose(ad_eng)
      self.removeIOobject(ad_ioNm)
    return var
  #[ .......... end of varReadSteps .............. ]#

  #[ Generate the Fourier space grid (note that for colorplots
  #[ matplotlib expects a nodal grid).
  def kGrid(self, **kwargs):
    Nekx  = self.attrRead('Nx', **kwargs)
    kxMin = self.attrRead('dx', **kwargs)
    dim   = np.size(Nekx)
    Nkx   = (Nekx+1)//2;  Nkx[1] = Nekx[1]; #[ Nekx[1] doesn't include negative ky's.
    kx    = [ np.array([kxMin[d]*i for i in range(Nekx[d])]) for d in range(dim) ]
    #[ Negative kx modes in increasing order:
    for i in range(Nkx[0],Nekx[0]):
      kx[0][i] = -(Nkx[0]-1-(i-Nkx[0]))*kxMin[0]

    if 'nodal' in kwargs:
      kxN = [ np.concatenate(([kx[d][0]-kxMin[d]/2.],
                              (0.5*(kx[d][:-1]+kx[d][1:])).tolist(),
                              [kx[d][-1]+kxMin[d]/2.])) for d in range(dim)]
      if kwargs['nodal']:
        return kxN
      else:
        return kx
    else:
      return kx
  #[ .......... end of kGrid .............. ]#

  #[ Generate the real space grid (note that for colorplots
  #[ matplotlib expects a nodal grid).
  def xGrid(self, **kwargs):
    fileGridType = self.attrRead('grid_type', **kwargs)
    if fileGridType == 'MUGY_REAL_GRID':
      Nx  = self.attrRead('Nx', **kwargs)
      dx  = self.attrRead('dx', **kwargs)
      dim = np.size(Nx)
    else:
      #[ Reconstruct real grid based on the Fourier grid in the file.
      Nekx  = self.attrRead('Nx', **kwargs)
      kxMin = self.attrRead('dx', **kwargs)
      Nkx = (Nekx+1)//2;  Nkx[1] = Nekx[1]; #[ Nekx[1] doesn't include negative ky's.
      Nx  = 2*(Nkx-1)+1;
      dim = np.size(Nx)
      Lx = [ 2.*np.pi/kxMin[d] for d in range(dim)]
      dx = [ Lx[d]/max(1.,Nx[d] - Nx[d] % 2) for d in range(dim)];

    x = [ np.array([(i-Nx[d]*0.5)*dx[0] for i in range(Nx[d])]) for d in range(dim) ]

    if 'nodal' in kwargs:
      xN = [ np.concatenate(([x[d][0]-dx[d]/2.],
                             (0.5*(x[d][:-1]+x[d][1:])).tolist(),
                             [x[d][-1]+dx[d]/2.])) for d in range(dim)]
      if kwargs['nodal']:
        return xN
      else:
        return x
    else:
      return x
  #[ .......... end of xGrid .............. ]#

  #[ Calculate the FLR operators.
  #[ 'tauIn'=Ti0/Ts0 and 'muIn'=sqrt(m_i/m_s), where
  #[ Ti0/m_i are the temperature and mass of the reference ions,
  #[ and Ts0/m_s are the temperature and mass of the species
  #[ whose FLR operators are desired.
  def calcFLR(self, kxIn, tauIn, muIn, **kwargs):
    kxSq    = [np.power(kxIn[0],2), np.power(kxIn[1],2)]
    kperpSq = np.zeros((np.size(kxIn[0]),np.size(kxIn[1])))
    for i in range(np.size(kperpSq,0)):
      for j in range(np.size(kperpSq,1)):
        kperpSq[i,j] = np.add(kxSq[0][i],kxSq[1][j])
    kperp = np.sqrt(kperpSq)
  
    krho = tauIn*kperp/muIn

    b   = np.power(krho,2)
    bSq = np.power(b,2)
  
    Gamma0 = scsp.ive(0,b)
    Gamma1 = scsp.ive(1,b)
    avgJ0  = np.sqrt(Gamma0)
  
    GammaRat  = Gamma1/Gamma0
    hatLap    = b*(GammaRat-1.0)
    hatLapSq  = np.power(hatLap,2)
    hathatLap = b*(0.5*GammaRat-1.0)-0.25*bSq*(3.0+GammaRat)*(GammaRat-1.0)
  
    Db = 1.0+hathatLap-0.25*hatLapSq
    Nb = 1.0+hathatLap-0.5*hatLapSq
  
    Sb = (avgJ0*Nb)/Db
  
    flrDict = {}
    if 'only' in kwargs:
      #[ Output only desired arrays (saves memory).
      for iSt in kwargs['only']:
        if iSt == "Gamma0":
          flrDict["Gamma0"]    = Gamma0
        elif iSt == "avgJ0":
          flrDict["avgJ0"]     = avgJ0
        elif iSt == "hatLap":
          flrDict["hatLap"]    = hatLap
        elif iSt == "hathatLap":
          flrDict["hathatLap"] = hathatLap
        elif iSt == "Db":
          flrDict["Db"]        = Db
        elif iSt == "Nb":
          flrDict["Nb"]        = Nb
        elif iSt == "Sb":
          flrDict["Sb"]        = Sb
    else:
      #[ Output all FLR functions.
      flrDict["Gamma0"]    = Gamma0
      flrDict["avgJ0"]     = avgJ0
      flrDict["hatLap"]    = hatLap
      flrDict["hathatLap"] = hathatLap
      flrDict["Db"]        = Db
      flrDict["Nb"]        = Nb
      flrDict["Sb"]        = Sb
    return flrDict
  #[ ........... End of calcFLR ................... ]#

  #[ Shift grid with the kx's going from 0 to the largest
  #[ positive kx followed by the negative kx's in
  #[ ascending order, to a space with kx=0 in the center.
  #[ Also return the number elements a variable needs to be
  #[ rolled by in order to center it at kx=0 as well.
  def kx0centerGrid(self, kxIn):
    kx = np.copy(kxIn[0])
    Nekx = np.size(kx[0])

    #[ Find max positive kx.
    Nkx = 0
    for i in range(Nekx):
      if kx[i] < 0.:
        break;
      else:
        Nkx += 1

    return np.append(kx[Nkx:], kx[:Nkx]), Nekx-Nkx 
  #[ ........... End of kx0centerGrid ................... ]#

  #[ Shift a variable defined on kx-ky space with the kx's
  #[ going from 0 to the largest positive kx followed by the
  #[ negative kx's in ascending order; to a space with kx=0
  #[ in the center.
  def kx0centerVar(self, varIn, **kwargs):
    if 'kx' in kwargs:
      _, nroll = self.kx0centerGrid(kwargs['kx'])
    elif 'rollnum' in kwargs:
      nroll = kwargs['rollnum']
    else:
      sys.exit("kx0centerVar: must provide either 'kx' or 'rollnum' arguments.")
    return np.roll(var, rollnum, axis=0)
  #[ ........... End of kx0centerVar ................... ]#
