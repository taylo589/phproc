'''
phmovie.py - a python module for processing frames of phase images and displaying them with a movie
'''
# plotting crap
#import matplotlib
#matplotlib.use('TkAgg')

import os
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from datetime import datetime
from phclasses import *



class Movie:
  def __init__(self,Dataset,iFrame=0,R=[-30,-20,-30,-20],CLIM=[-3.,3.]):
    self.REFREGION = (slice(R[0],R[1]),slice(R[2],R[3]))
    self.CLIM = CLIM
    self.retrieval_method = Dataset.retrieval_method
    self.FrameList = Dataset.FrameList
    rFrame = Frame(Dataset.FrameList[iFrame])
    self.phim0,self.arg1,self.arg2 = self.retrieval_method(rFrame.data,ARGS=True)
#    self.clim = self.implt.get_clim()


  def animation_function(self,FRAMEFILE):
    iFrame = Frame(FRAMEFILE)
    phim = self.retrieval_method(iFrame.data,self.arg1,self.arg2)
    phim_i0 = phim - self.phim0
    phim_i0r,bg = phplot.sub_ref(phim_i0,self.REFREGION)
    self.bg = self.a * bg + (1. - self.a)*self.bg

    self.implt.set_data(phim_i0r)
    self.fig.suptitle(FRAMEFILE[-8:-4])



  def play(self):
    self.fig = plt.figure()
    self.n = 1.
    self.a = 1./(5.+1.)
    phim0r,self.bg = phplot.sub_ref(self.phim0,self.REFREGION)
    self.implt = phplot.imageplot(phim0r,datalim=self.CLIM)
    ani = anim.FuncAnimation(self.fig,self.animation_function,self.FrameList,\
        interval=300,repeat=False,blit=False)
    plt.show()

# animation utility function.
#def phaseani(FILENAME,FILTER,REFREGION,bias,TIMESTAMP):
#  global implot
#  global ax
#
#  # read interferogram file
#  im,timestamp = util.im_read(FILENAME,0.0,TIMESTAMP)
#  # process interferogram and get phase
#  im = methods.get_phase(im,gMETHOD,FILTER,REFFIELD=field,LEVEL=False)
##  im = methods.hilbert_phase(im,pre_filter=FILTER,LEVEL=False)
##  im = methods.arccos_phase(im,pre_filter=FILTER,LEVEL=False)
#  # subtract bias plane taken from inital phase frame
#  im -= bias
#  # subtract from reference region of frame - that is - make
#  # a relative phase measurement to the frame itself
##  im -= np.mean(np.mean(im[-10:-1,0:10]))
##  im[10:20,10:20] = 10.
####
#  im_REF = np.mean(np.mean(im[REFREGION]))
#  im -= im_REF
#  im[REFREGION] = 10.
####
##  im = ndi.gaussian_filter(im,5)
#
#  # add timestamp and update plot
##  im = util.add_timestamp(im,timestamp)
#  implot.set_data(im.T)
#  ax.set_title(util.get_name(FILENAME) + '\n' + str(im_REF))
#
#def ph_movie(imfiles,METHOD,FILTER,RANGE,REFREGION,REFFIELD,SAVE,PATH,TIMESTAMP=True):
#  global implot
#  global ax
#  global gMETHOD
#  global field
#   
## read initial fame
#  im0,timestamp = util.im_read(imfiles[0],0.0,TIMESTAMP)
#
## see if image recognition works
##  normts = timestamp + np.min(np.min(timestamp))
##  normts += 100
##  normts **= 1.2
##  normts = normts/np.max(np.max(normts))
##  pilts = fromarray(np.uint8(normts*255))
##  plt.imshow(pilts)
##  plt.colorbar()
##  plt.show()
##  print(pytesseract.image_to_string(pilts))
##  return
#
## process initial frame, saving level-plane for future leveling
#  gMETHOD = METHOD # set global method
#  if REFFIELD is not None:
#    field = ndi.imread(REFFIELD)
#
#  im0_u = methods.get_phase(im0,METHOD,FILTER,REFFIELD=field,LEVEL=False)
##  im0_u = methods.hilbert_phase(im0,pre_filter=FILTER,LEVEL=False)
##  im0_u = methods.arccos_phase(im0,pre_filter=FILTER,LEVEL=False)
#  im0_l = level.diff_level(im0_u)
#  im0_ph = im0_u - im0_l
#  im0_REF = np.mean(np.mean(im0_ph[REFREGION]))
#  im0_ph -= im0_REF
#  im0_bias = im0_l #+ im0_ph # add this initial frame to the bias to be subtracted
#  #im0_bias += ndi.gaussian_filter(im0_ph,5)
####
#  im0_bias += im0_ph
####
#  im0_ph = ndi.gaussian_filter(im0_ph,1)
## prepare plot
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
##  implot = plt.imshow(im0_ph,cmap=plt.cm.gray,extent=[0,1,0,0.5])
#  #yext = float(im0_ph.shape[0])/im0_ph.shape[1]
#  yext = float(im0_ph.shape[1])/im0_ph.shape[0]
#  #implot = plt.imshow(im0_ph,cmap=plt.cm.gray,extent=[0,1,0,yext])
#  implot = plt.imshow(im0_ph.T,cmap=plt.cm.gray,extent=[0,1,0,yext])
##  implot = plt.imshow(im0_ph,cmap=plt.cm.hot,extent=[0,1,0,yext])
#  DIM = settings.D_CAM*1e6
#  plt.xlabel("x / {0} um".format(DIM))
#  plt.ylabel("y / {0} um".format(DIM))
#  plt.colorbar()
#  implot.set_clim(RANGE)
#
## call animation function
#  ani = anim.FuncAnimation(fig,phaseani,imfiles,\
#  fargs=(FILTER,REFREGION,im0_bias,implot),interval=300,repeat=False,blit=False)
#
#  if SAVE:
#    FILESTAMP = datetime.today().strftime("%Y%m%d%H%M")
#    if os.path.exists('/usr/local/bin/ffmpeg'):
#      ani.save(PATH + 'phase_' + FILESTAMP + '.mp4',writer='ffmpeg',bitrate=1000)
#    else:
#      ani.save(PATH + 'phase_' + FILESTAMP + '.mp4',writer='avconv',fps=5,bitrate=10000)
#  else:
#    plt.show()
