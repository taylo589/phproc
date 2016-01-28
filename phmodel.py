'''
phmodel.py - physical model for sample heat transfer

Revisions:
25/1/2016 - LNT - first
'''

import numpy as np
from scipy.special import erf
import phplot

class temperature_field(object):

  # impulse response from Carslaw & Jaeger
  # centered at "x", "y", "z"
  def T_impulse(self,A,a,x,y,z,tau):
    return A*np.exp(-((self.X-x)**2+(self.Y-y)**2+(self.Z-z)**2)/(4*a*(tau)))

  # symbolic model transferred from Mathematica
  # centered on "x"
  def T_symb(self,A,a,x,tau):
    X = self.X[:,:,0]
    return A*np.pi*\
        (2*np.exp(-(X-x)**2/(4*a*tau))-\
        2*np.sqrt(np.pi)*(np.abs(X-x)-\
        (X-x)*\
        erf(X/np.sqrt(tau))))

  def __init__(self):
    nx,ny,nz = 100,100,100
    lx,ly,lz = 1.,1.,1.
    x,y,z = np.linspace(-lx,lx,nx),np.linspace(-ly,ly,ny),np.linspace(-lz,0,nz)
    self.X,self.Y,self.Z = np.meshgrid(x,y,z)
    nt = 100
    lt = 10
    self.t = np.linspace(0.,lt,nt)
    self.I = self.T_impulse(1.,1.,0.,0.,0.,self.t[1])
    phplot.imageshow(np.sum(self.I,2))
    self.S = self.T_symb(1.,1.,0.,self.t[99])
    phplot.imageshow(self.S)
