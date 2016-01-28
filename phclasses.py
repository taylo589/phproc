'''
phclasses.py - classes useful to the methods used to process frames

mainly consists of:
  Frameset - the data set of frames grabbed from the camera in the setup
  Frame - a single-frame object. ideally contains everything one would need when performing operations on a frame
'''
# useful imports
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
import glob
import re
import json
import phretrieve
import phplot

# class comprising data and methods for the data set as a whole
class Frameset(object):

  # useful strings for the methods - not necessesary for user access
  DATAERROR_STR = "Incomplete experimental parameter data. Some functionality may not be present. Consider using a different data set."

  # arguments:
  # FRAMEPATH - path to the folder where the camera frames are stored
  # FILEEXT - file extension string for the camera frame files
  # PROCMETHOD - method to use to retrieve the phase of the image
  # SAMPLEPATH - path to a frame showing only the sample image (no interference fringes)
  def __init__(self,FRAMEPATH=None,FRAMESIZE=1e-3,FILEEXT ='.tif',retrieval_method=None,SAMPLEPATH=None):

    # store the quantities useful to the frameset
    self.FRAMEPATH = FRAMEPATH
    self.FrameList = self.sort_files(FRAMEPATH,FILEEXT)
    self.retrieval_method = retrieval_method

    # get the sample-arm-only frame as a Frame object
    if SAMPLEPATH is not None:
      self.SAMPLEPATH = SAMPLEPATH
      self.SampleFrame = Frame(SAMPLEPATH)

    # try to get the experimental parameters
    try: 

      self.exp_params = self.JSONload(FRAMEPATH + "exp_params.json")
      self.camsettings = self.JSONload(FRAMEPATH + "camsettings.json")
      self.rdata = self.JSONload(FRAMEPATH + "data.json")
    except IOError as ioe:
      print ioe
      print self.DATAERROR_STR
    except ValueError as ve:
      print ve
      print self.DATAERROR_STR
    finally:
      self.dimension = FRAMESIZE # meters (approximate - check calibration frame)

  # load from JSON file from path
  def JSONload(self,PATH):
    with open(PATH,'r+') as f:
      return json.load(f)

  # intelligently gather all the frame files on path
  def sort_files(self,PATH,FILEEXT):
    FILENAMES = glob.glob(PATH + '*' + FILEEXT)
    FILENAMES = sorted(FILENAMES,\
        key=lambda filename: [int(c) for c in re.split('(-*\d+\.\d*)',filename) if c.isdigit()])
    return FILENAMES



# class comprising data and methods for an individual frame in the data set
class Frame:
  def __init__(self,FRAMEFILE=None):
    self.FRAMEFILE = FRAMEFILE
    self.data = self.get_data(self.FRAMEFILE)
    self.retrieve_default = phretrieve.hilbert_retrieval

  def get_data(self,FRAMEFILE):
    framedata = ndi.imread(FRAMEFILE)
    return np.array(framedata,dtype=np.float_)

  def set_phase(self,retrieve_function=None,REFREGION=None,ARGS=False):
    if retrieve_function is not None:
      self.phase = retrieve_function(self.data,ARGS=ARGS)
    else:
      self.phase = self.retrieve_default(self.data,ARGS=ARGS)

  def min2D(self,data):
    return np.min(np.min(data))

  def max2D(self,data):
    return np.max(np.max(data))

  def get_ext(self):
    return [0,np.shape(self.data)[0],0,np.shape(self.data)[1]]

  def get_lim(self):
    return min2D(self.data),max2D(self.data)

  def showframe(self):
    phplot(self.data)
    plt.show()

