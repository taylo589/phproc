'''
phplot.py - a collection of useful plotting functions for plotting image data

'''
import matplotlib.pyplot as plt
import numpy as np

def imageplot(data,dataextent=None,datalim=None):
  if dataextent is None:
    dataextent = [0,1,0,1]
  if datalim is None:
    iplt = plt.imshow(data,cmap=plt.cm.gray,extent=dataextent)
  else:
    iplt = plt.imshow(data,cmap=plt.cm.gray,extent=dataextent,clim=datalim)
  plt.colorbar()
  return iplt

def imageshow(data,dataextent=None,datalim=None):
  imageplot(data,dataextent,datalim)
  plt.show()

def dataplot(*data):
  for dataset in data:
    plt.plot(dataset)

def datashow(*data):
  dataplot(*data)
  plt.show()

def dB_data(data):
  return np.log1p(np.abs(np.fft.fftshift(data)))

def spectrogram(data):
  return dB_data(np.fft.fft2(data))

def sub_ref(phase,REFREGION,MARKER = 1e3):
  ph = np.copy(phase)
  bg = np.mean(np.mean(ph[REFREGION]))
  ph -= bg
  ph[REFREGION] = MARKER
  return ph,bg

def plot_hilbert_steps(Dataframe,SAVE=False):
  data,f_data,ff_data,h_data,w_data,u_data,l_data,l_plane,ph_data\
   = phretrieve.hilbert_phase(Dataframe.data,STEPS=True)
  steps = (data,\
           self.dB_data(f_data),\
           self.dB_data(ff_data),\
           self.dB_data(h_data),\
           w_data,\
           u_data,\
           l_data,\
           l_plane,\
           ph_data)
  for i,step in enumerate(steps):
    plt.figure(i)
    steplim = (Dataframe.min2D(step),Dataframe.max2D(step))
    imageplot(step,datalim=steplim)
    plt.colorbar()

  if not SAVE:
    plt.show()

  else:
    print "Not implemented yet"
