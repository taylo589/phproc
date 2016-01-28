'''
phretrieve.py - a collection of methods used for phase processing and retrieval.

21/12/2015 - LNT - first.
'''
# useful imports
import numpy as np
import scipy.ndimage as ndi
from scipy.signal import correlate2d
from skimage.restoration import unwrap_phase
# for testing purposes
import phplot

# subtract a reference field from an interferogram
# (helps eliminate non-phase background noise)
def subtract_field(data,ref_data):
  # subtract the max point of the data
  data_max = np.max(np.max(data))
  data_n = data - data_max

  phplot.imageshow(data_n)
  phplot.imageshow(phplot.spectrogram(data_n))
  # subtract the reference field
#  data_n = data_n - ref_data
  wx = np.hanning(data.shape[0])
  wy = np.hanning(data.shape[1])
  w2 = np.sqrt(np.outer(wx,wy))
  phplot.imageshow(data_n)
  phplot.imageshow(phplot.spectrogram(data_n*w2))

  # set 0-points of reference data to some non-0 value (median)
  # this is important for the division step
  ref_data[ref_data == 0] = np.median(ref_data)

  # divide by reference field
  data_n /= 2*np.sqrt(ref_data*data_max)

  return data_n

# partition field to remove 0-sections (they confuse the phase algorithm)
# this is frame specific and designed for 2 metal lines blocking the field
# as was the case for measurements made in Fall 2015 with DC Chien-type 
# phase-imaging samples
def partition_field(data,ref_data):

  # get the derivative of data (Sobel-filter)
  sobel_data = ndi.sobel(ref_data)

  # sum along y-dimension (reduce the dimensionality of data)
  sum_sobel = np.sum(np.abs(sobel_data),0)
  # find the mean
  mean_sobel = np.mean(sum_sobel)
  # find the standard deviation
  std_sobel = np.std(sum_sobel)

  # find the max-points of Sobel-field
  # - they indicate transitions to zero field
  sum_sobel_n1 = sum_sobel[:-2] # negative 1 place
  sum_sobel_m = sum_sobel[1:-1] # middle place
  sum_sobel_p1 = sum_sobel[2:]  # positive 1 place
  threshold = 20e3

  # return the transition indecies (might need to add 1)
  indecies = np.where(sum_sobel_m > sum_sobel_n1 and\
      sum_sobel_m > sum_sobel_p1 and\
      sum_sobel_m > threshold)

  # check to see if we have 4 Sobel points where the field-dark transitions occur
  # this is frame specific
  if len(indecies) == 4:
    indecies[0] -= 3
    indecies[1] += 2
    indecies[2] -= 1
    indecies[3] += 2

    # grab the "clear data" should have no dark regions
    data_clear = data[:,indecies[3]:]
    # stitch the data around the dark regions
    data_stitched = np.concatenate((\
        data[:,:indecies[0]],\
        data[:,indecies[1]:indecies[2]],\
        data[:,indecies[3]:]),1)
    return data_clear,data_stitched
  else:
    print "Partitioning didn't work. Returning original field."
    return data

def min2D(data):
  return np.min(np.min(data))

def max2D(data):
  return np.max(np.max(data))

def get_band(data,wave=None,bw=0.5):
  if wave is None:
    wave = get_wave(data)
  nx,ny = np.shape(data)
  kmax = np.max(wave)
  imax = np.where(wave == kmax)[0]
  nmax = (nx,ny)[imax]
  kadj = float(nmax/2. - kmax)/(nmax/2.)
#  print kadj
  return [(1.-bw)*kadj,(1.+bw)*kadj]


# peak-finding algorithm - not robust
# not well-tested
# a work-in-progress
# relies on DC term at center (FFT shift)
def get_disparate_peaks(data,a=1.):
  datamax = max2D(data)
  x_max,y_max = np.where(data == datamax)
  # minimum distance threshold for peak
  nx,ny = np.shape(data)
  n = np.sqrt(nx*ny)
  m = n*0.3

  while True:
    thresh = datamax*a
    x_peaks,y_peaks = np.where(data > thresh)
    d = (x_peaks - x_max)**2 + (y_peaks - y_max)**2
    n = len(d[d > m])
#    print d

    # evaluate peaks
    # not enough peaks, decrease threshold
    if n < 1:
      a *= 0.9
    # missed the peak - too many points, add to threshold factor
    elif n > 2:
      a *= 1.2
    # case where there is 1 unique peak outside of DC (mirrored in negative space)
    else:
#      print a 
#      print d
#      print x_peaks,y_peaks
      return x_peaks[-1],y_peaks[-1]
  
# find the k-vector of interferogram data
# first do the hilbert transform to get the phase!
def get_wave(data,DOFFT=False):
  if DOFFT:
    f_data = np.fft.fft2(data)
  else:
    f_data = np.copy(data)
  # find 2nd gradients
  gf_data = np.gradient(np.abs(np.fft.fftshift(f_data)))
  #gf_data = np.gradient(np.abs(f_data))
  sgf_data = +gf_data[0] + gf_data[1]
  gsgf_data = np.gradient(sgf_data)
  sgsgf_data =-( gsgf_data[0] + gsgf_data[1])
  
  #print np.where(ggf_data < 0)
#  phplot.imageshow(sgsgf_data)
  std = np.std(sgsgf_data)
  mean = np.mean(sgsgf_data)
  mx = max2D(sgsgf_data)
  a = 5.  # currently find peaks based on a threshold of max/5 
  # not a great way, because all data will have different values ....
  Kx,Ky = get_disparate_peaks(sgsgf_data)
  nx,ny = np.shape(data)
  #Kx,Ky = np.where(sgsgf_data > mx/a)
  K = (nx+(nx/2-Kx),ny/2 - Ky)
#  K = Kx,Ky
  f_data[K[0],K[1]] = 1e10
#  K = (K[0] + nx/2.), (K[1] + ny/2.)
#  Kabs = np.sqrt(K[0]**2 + K[1]**2)
#  k = K[0]/Kabs,K[1]/Kabs

  
#  phplot.imageshow(data)
#  phplot.imageshow(phplot.dB_data(f_data))
  return K



# a leveling function copied directly from phaseproc_0_2 level.py
def diff_level(imdata,K=False):
# assuming data is dominated by a level plane, find derivatives and 
# extrapolate them to a plane
  ky,kx = np.gradient(imdata)
  im00 = imdata[0,0]
#  print np.mean(kx),np.mean(ky)
  kx = np.mean(kx)
  ky = np.mean(ky)
  if K:
    return kx,ky,im00
  else:
    nx = np.size(imdata,1)
    ny = np.size(imdata,0)
    x,y = np.meshgrid(np.linspace(0,nx-1,nx),np.linspace(0,ny-1,ny),
    sparse=False, indexing='xy')
    return kx*x+ky*y + im00

# simple bandpass fiter in fourier-domain
def f_bandpass(f_data,band):
  filt_data = np.copy(f_data)

  # perform the bandpass as 2 lowpass operations
  filt_data = ndi.fourier_gaussian(filt_data - ndi.fourier_gaussian(filt_data,band[1]),band[0])
  #phplot.imageshow(phplot.dB_data(filt_data))

  return filt_data

# "hilbert transform" function (in fourier-space)
def hilbert_transform(f_data,direction='x'):
  # do a copy of the data (just in case)
  h_data = np.copy(f_data)

  # multiply image field with step function in x or y
  if direction == 'x':
    x0 = np.around(np.size(h_data,1)/2)
    h_data[:,x0:] *= 1.
    h_data[:,:x0] *= 0.
  elif direction == 'y':
    y0 = np.around(np.size(h_data,0)/2)
    h_data[y0:,:] *= 0.
    h_data[:y0,:] *= 1.

  # return the filtered spectral data
  return h_data

def f_gaussian(data,k,s=0.05):
  nx,ny = np.shape(data)
  N = max([nx,ny])
  KY,KX = np.meshgrid(np.arange(nx),np.arange(ny),sparse=False,indexing='xy')
  K = (nx-k[0])+nx/2,ny/2-k[1]
#  print K
#  data = np.fft.fftshift(data)
  #g = 1e6*np.exp(-((KX-K[0])**2+(KY-K[1])**2)/(2*(s*N)**2))
  #return np.fft.ifftshift(data+g)
  g = np.exp(-((KX-K[0])**2+(KY-K[1])**2)/(2*(s*N/2)**2))
#  phplot.imageshow(phplot.dB_data(np.fft.ifftshift(g)))
#  phplot.imageshow(np.abs(np.fft.ifft2(np.fft.ifftshift(g))))
  return np.fft.ifftshift(g)


# produces a wave with X and Y mesh data
def wave_function(X,Y,K,phi=0.):
  return np.cos(K[0]*X*2*np.pi + K[1]*Y*2*np.pi + phi)

def sum2D(data):
  return np.sum(np.sum(data))

def mean2D(data):
  return np.mean(np.mean(data))

def correlation_retrieval(data,K=None,ARGS=False):
  f_data = np.fft.fft2(data)
  kx = np.fft.fftfreq(f_data.shape[0])
  ky = np.fft.fftfreq(f_data.shape[1])
  # get k-vector, if not give
  if K is None:
    K = get_wave(f_data)

  # get the frequency data information
  Kfreq = kx[K[0]],ky[K[1]]
  # make an X-Y grid
  pi = np.pi
  nx,ny = np.shape(data)
  Y,X = np.meshgrid(np.arange(ny),np.arange(nx),\
      sparse=False,indexing='xy')
  print Kfreq
  # produce a correlation function
  phplot.imageshow(wave_function(X,Y,Kfreq))
  zsum = np.zeros(data.shape)
  Nphi = 1.
  a_phi = np.linspace(0,2*np.pi,Nphi)
  for phi in a_phi:
    zsum += correlate2d(data,wave_function(X,Y,K,phi),mode='same')

  phplot.imageshow(zsum)
  return data
    
def wave_retrieval(data,r_data=None,filt=None,STEPS=False,s=0.1,ARGS=False):
  # go to fourier space
  f_data = np.fft.fft2(data)
  # get the k-vector (lots of finagling here for fft shift crap)
  if r_data is None:
    K = get_wave(f_data)
    kx,ky = np.fft.fftfreq(f_data.shape[0]),np.fft.fftfreq(f_data.shape[1])
    k = kx[K[0]],ky[K[1]]
    nx,ny = data.shape
    Y,X = np.meshgrid(np.arange(nx),np.arange(ny),sparse=False,indexing='xy')
  #  phplot.imageshow(data)
  #  phplot.imageshow(np.real(r_data))
  #  phplot.imageshow(np.real(r_data*np.exp(+1j*2.*np.pi*(k[0]*X+k[1]*Y))))
    # reference beam
    r_data = np.exp(+1j*2.*np.pi*(k[0]*X+k[1]*Y))


  # filter k-space with a gaussian around data of interest
  if filt is None:
    try:
      gf_data = f_gaussian(f_data,K)*f_data
    except UnboundLocalError:
      K = get_wave(f_data)
    finally:
      gf_data = f_gaussian(f_data,K)*f_data
      filt = np.ones(f_data.shape)
      filt = f_gaussian(filt,K)*filt
#      band = get_band(data,K)
#      filt = f_bandpass(filt,band)
#      filt = hilbert_transform(filt,'x')
  else:
    gf_data = filt*f_data
  # check it out!!
#  phplot.imageshow(phplot.dB_data(gauss_filter))
  # go back to xy space
  g_data = np.fft.ifft2(gf_data)

  o_data = g_data*r_data  # retrieve the wrapped phase from the digital reference beam
  wrapped_phase = np.arctan2(o_data.imag,o_data.real)
  # unwrap the phase (may be unneccesary)
  unwrapped_phase = unwrap_phase(wrapped_phase)
  phase = unwrapped_phase
  #phase = wrapped_phase

  if STEPS:
    return data,f_data,gauss_filter,wrapped_phase,unwrapped_phase,plane,phase

  if ARGS:
    return -phase,r_data,filt

  return -phase

# main hilbert phase-retrieval function
# 2nd argument is generic 2ndary data
def hilbert_retrieval(data,plane=None,filt=None,STEPS=False,direction='x',ARGS=False):
  # fourier transform data and filter
  window_f = np.hanning
  w2 = np.sqrt(np.outer(window_f(data.shape[0]),window_f(data.shape[1])))
  f_data = np.fft.fft2(data)
  if filt is None:
    band = get_band(f_data)
    filt_f_data = f_bandpass(f_data,band)
    hilbert_f_data = hilbert_transform(filt_f_data,direction)
    filt = np.ones(f_data.shape)
    filt = f_bandpass(filt,band)
    filt = hilbert_transform(filt,direction)
  else:
    hilbert_f_data = f_data *filt

  # hilbert-transform
#  phplot.imageshow(phplot.dB_data(filt_f_data))
  hilbert_data = np.fft.ifft2(hilbert_f_data)

  # retrieve the wrapped phase
  wrapped_phase = np.arctan2(hilbert_data.imag,hilbert_data.real)

  # unwrap the phase
  unwrapped_phase = unwrap_phase(wrapped_phase)

  # flatten the data (if required - some methods can use the same level for all data)
  if plane is None:
    plane = diff_level(unwrapped_phase)
  phase = unwrapped_phase - plane

  if STEPS:
    return data,f_data,filt_f_data,hilbert_f_data,hilbert_data,wrapped_phase,unwrapped_phase,plane,phase

  if ARGS:
    return -phase,plane,filt
  return -phase

def plane_wave(X,Y,kx,ky):
  return np.exp(-1j*2*np.pi*((X*kx)+(Y*ky)))
  
def gaussian(X,Y,mu_x,mu_y,sigma_x,sigma_y=None):
  if sigma_y is None:
    sigma_y = sigma_x
  return np.exp(-((X-mu_x)**2)/(2*sigma_x**2)-((Y-mu_y)**2)/(2*sigma_y**2))
  
def find_carrier(data,n=None):
  if n is None:
    nx,ny = data.shape[0]-1,data.shape[1]
#  window_f = np.hamming
#  w2 = np.sqrt(np.outer(window_f(data.shape[0]),window_f(data.shape[1])))
  f_data = np.fft.rfft2(data)
  Kx,Ky = np.fft.rfftfreq(nx),np.fft.fftfreq(ny)
  print f_data.shape
  ik = get_disparate_peaks(f_data)
#  print ik
#  f_data[ik] = 1e10
#  phplot.dBshow(np.fft.fftshift(f_data,0),SHIFT=False)
  kx,ky = Kx[ik[1]],Ky[ik[0]]
#  print kx,ky
  return kx,ky

# find a gaussian centered on the wave
def find_gaussian(data,kx,ky):
  # check out a gaussian in the area
  Kx,Ky = np.fft.rfftfreq(data.shape[0]-1),np.fft.fftfreq(data.shape[1])
  sigma = 20. * np.max(Ky)/256
  KX,KY = np.meshgrid(Kx,Ky)
  G0 = gaussian(KX,KY,kx,ky,sigma)
  # combine gaussian filter with the rest of the space
  G = np.fft.ifftshift(\
    np.concatenate(\
    (np.zeros(G0.shape),\
    np.fft.fftshift(gaussian(KX,KY,kx,ky,sigma),0)),1))
  return G

def wave_retrieval_2(data,ARGS):
  kx,ky = find_carrier(data)
  G = find_gaussian(data,kx,ky)
  I = np.fft.fft2(data)*G
#  phplot.dBshow(I)
  nx,ny = data.shape
  X,Y = np.meshgrid(np.arange(nx),np.arange(ny))
  o = np.fft.ifft2(I)*plane_wave(X,Y,kx,ky)
  return np.arctan2(o.imag,o.real)

def compressed_wave_retrieval(data,ARGS):
  window_f = np.hamming

  kx,ky = find_carrier(data)
  nx,ny = data.shape
  X,Y = np.meshgrid(np.arange(nx),np.arange(ny))

  G = find_gaussian(data,kx,ky)
#  print float(ik[1])/256
#  print np.around(float(ik[1])/ax) - float(ik[1])/ax
  a = [11./256.,12./256.]
  G_kernel = ndi.zoom(G,a)
#  phplot.dBshow(G_kernel)
  g_kernel = np.fft.ifft2(G_kernel)
  wzoom = np.sqrt(np.outer(window_f(g_kernel.shape[0]),window_f(g_kernel.shape[1])))
  g_kernel *= wzoom
#  print g_kernel.shape
#  phplot.imageshow(g_kernel.real)
#  nx,ny = data.shape
#  X,Y = np.meshgrid(np.arange(nx),np.arange(ny))
  i = (ndi.convolve(data,g_kernel.real) + 1j*ndi.convolve(data,g_kernel.imag))
#  phplot.dBshow(np.fft.fft2(i))
  r = plane_wave(X,Y,kx,ky)
  o = i*r
  ph = np.arctan2(o.imag,o.real)
#  phplot.imageshow(p)
  return ph

# another method
def holographic_retrieval(data,ARGS):
  # find a window for the fft
  window_f = np.hamming
  w2 = np.sqrt(np.outer(window_f(data.shape[0]),window_f(data.shape[1])))
  # get the wave-vector
  kx,ky = find_carrier(data)
  Kx,Ky = np.fft.fftfreq(data.shape[0]),np.fft.fftfreq(data.shape[1])
  KX,KY = np.meshgrid(Kx,Ky)
  s2 = np.sqrt(kx**2 + ky**2)/4.
  G2 = gaussian(KX,KY,0,0,s2)
#  phplot.dBshow(G2)

  nx,ny = data.shape
  X,Y = np.meshgrid(np.arange(nx),np.arange(ny))
  r = plane_wave(X,Y,kx,ky)
  i2 = data*r
  I2 = np.fft.fft2(w2*i2)
  O2 = I2 * G2
  o2 = np.fft.ifft2(O2)
  p2 = np.arctan2(o2.imag,o2.real)
#  phplot.dBshow(I2)
#  phplot.dBshow(O2)
#  phplot.imageshow(p2)
  return p2

# convolution form of the method
def compressed_holographic_retrieval(data,ARGS):
  window_f = np.hamming
  kx,ky = find_carrier(data)
  Kx,Ky = np.fft.fftfreq(data.shape[0]),np.fft.fftfreq(data.shape[1])
  KX,KY = np.meshgrid(Kx,Ky)
  s2 = np.sqrt(kx**2 + ky**2)/4.
  G2 = gaussian(KX,KY,0,0,s2)

  nx,ny = data.shape
  X,Y = np.meshgrid(np.arange(nx),np.arange(ny))
  r = plane_wave(X,Y,kx,ky)
  i3 = data*r

#  print np.around(float(s2*256)/ax) - float(s2*256)/ax
  a3 = 14./256
  G3_kernel = ndi.zoom(G2,a3)
#  phplot.dBshow(G3_kernel)
  g3_kernel = np.fft.ifft2(G3_kernel)
  w3zoom = np.sqrt(np.outer(window_f(g3_kernel.shape[0]),window_f(g3_kernel.shape[1])))
  g3_kernel *= w3zoom
  o3 = (ndi.convolve(i3.real,g3_kernel.real) -\
        ndi.convolve(i3.imag,g3_kernel.imag) +\
        1j*ndi.convolve(i3.real,g3_kernel.imag) +\
        1j*ndi.convolve(i3.imag,g3_kernel.real))
  p3 = np.arctan2(o3.imag,o3.real)
#  phplot.imageshow(p3)
  return p3

# generic testing method for phase-retrieval
def test_retrieval(data,r=None,s=None,ARGS=False):
  if ARGS or r is None or s is None:
    kx,ky = find_carrier(data)
    nx,ny = data.shape
    X,Y = np.meshgrid(np.arange(nx),np.arange(ny))
    r = plane_wave(X,Y,kx,ky)
    s = 1./(1.8*np.sqrt((kx)**2+(ky)**2))
  i3 = data*r
  o3 = ndi.gaussian_filter(i3.real,s) + 1j*ndi.gaussian_filter(i3.imag,s)
  #phplot.dBshow(np.fft.fft2(ndi.gaussian_filter(data,s)))
  p3 = np.arctan2(o3.imag,o3.real)
  if not ARGS:
#    return unwrap_phase(p3)
    return p3
  else:
#    return unwrap_phase(p3),r,s
    return p3,r,s
