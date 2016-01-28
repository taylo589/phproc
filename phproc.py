'''
phproc.py - this file, as part of the phaseproc3 suite, is my first attempt at an object-oriented approach to some of the methods in phaseproc_0_1 and phaseproc_0_2. There were a lot of redundances in addition to the useful extra features in those modules. I wanted to clean things up and take a more object oriented approach.

I would like an object representing the entire data space - but this would be quite large: 256x256x400. I'm trying to figure out how to do this.

Conventions:
ALLCAPS - a string
alllowercase - a single variable, list, tuple or dictionary
Capitalized - a class instance

Revisions:
16/12/2015 - LNT - first
'''
import os
import scipy.ndimage as ndi
import json
import glob
# phaseproc3 imports
from phclasses import *
import phretrieve
import phmovie
import phmodel


def walk_path(PATH):
  FOLDERS = [f_name[0] +'/' for f_name in os.walk(PATH)]
  i = 0
  for FOLDER in FOLDERS:
    print "[{0}] \"{1}\"".format(i,FOLDER.split('/')[-2])
    i += 1
  selection = raw_input(\
  "Select a number. \'q\' to exit: ")
  try:
    return FOLDERS[int(selection)]
  except ValueError as ve:
    print ve
    if selection != 'q':
      print "Oops! Try again!"
      walk_path(PATH)
    else:
      print "You selected \'q\'. Quitting..."
      return selection

def ph_retrieve(Dataset):
  sFrames = [0,len(Dataset.FrameList)]
#  iframe = int(float(\
#      raw_input("Pick a frame from {0} to {1}: ".format(sFrames[0],sFrames[1]))))
  iframe = int(float('0'))
  rFrame = Frame(Dataset.FrameList[iframe])
#  data_n = phretrieve.subtract_field(rFrame.data,Dataset.SampleFrame.data)
  rFrame.set_phase(Dataset.retrieval_method)
  phplt = phplot.imageshow(rFrame.phase)
#  phdata2 = rFrame.get_phase(PH_RETRIEVAL['h'])
#  phplt = phplot.imageshow(phdata2)
#  tfield = phmodel.temperature_field()
  return True

def ph_steps(Dataset):
  print "Nothing here yet"
  return True

def ph_movie(Dataset):
  movie = phmovie.Movie(Dataset)
  movie.CLIM = Defaults.clim
  movie.REFREGION = Defaults.refregion
  movie.play()
  return True

def ph_temporal(Dataset):
  REFREGION = Defaults.refregion
  INTEGRATE = Defaults.integrate
  ref = np.zeros(len(Dataset.FrameList))
  rel = np.copy(ref)

  # loop through the data frames
  for i,framepath in enumerate(Dataset.FrameList):
    iFrame = Frame(framepath)
    phase = Dataset.retrieval_method(iFrame.data)
    if i == 0:
      phase0 = np.copy(phase)

    # subtract the first, null field
    phase -= phase0

    # track the reference region over time (without smoothing data)
    ref[i] = np.mean(np.mean(phase[REFREGION]))
    rel[i] = np.mean(np.mean(phase[INTEGRATE])) - ref[i]
#  t = np.arange(ref.shape)/Dataset.camsettings["FrameRate"]
  phase[REFREGION] = 1e4
  phase[INTEGRATE] = 1e4
  phplot.imageshow(phase,datalim=[-3.,3.])
  phplot.datashow(ref,rel,ref+rel)


def ph_quit():
  print "Quitting..."
  return False


# point the main_loop at the directory with all the datasets
# give it a dict of program parameters
def main_loop(TESTPATH,SAMPLEPATH=None):
  # loop through program options
  QUIT = False
  CONT = False
  while not QUIT:

    # see if we want to loop through data set again
    if not CONT:
      PATH = walk_path(TESTPATH) # path returns 'q' if there is a problem
    # see if walk_path produced a quitting request
    if PATH == 'q':
      QUIT = True
      break
    print "You've chosen the data set at {0}".format(PATH)
    print "Operations available:"

    # print out operations and descriptions
    for key in OPERATIONS:
      print "[{0}] - {1}".format(key,OPERATION_DESC[key])

    # all operations should return True except the quit operation (confusing I know)
    OP = raw_input("Selection [q]: ")
    if OP not in OPERATIONS.keys():
      OP = 'q'

    print "Retrieval methods available:"
    # get the phase-processing method
    for key in PH_RETRIEVAL:
      print "[{0}] - {1}".format(key,PH_RETRIEVAL_DESC[key])
    # get the selection (defaults to hilbert retrieval)
    PH = raw_input("Selection [h]: ")
    if PH not in PH_RETRIEVAL.keys():
      PH = 'h'

    # load the Frameset object
    Dataset = Frameset(FRAMEPATH=PATH,SAMPLEPATH=SAMPLEPATH,retrieval_method=PH_RETRIEVAL[PH])

    # load the data set into the operation (maybe not all operations need this tho...)
    QUIT = not OPERATIONS[OP](Dataset)
    if QUIT:
      break

    # figure out what to do after we are finished
#    cont = raw_input("Operation complete. Hit \'q\' to quit, \'r\' to start over or ENTER to do another operation: ")
    # TEST input
    cont = 'q'
    if cont == 'q':
      QUIT = True
      break
    CONT = cont !='r'

# phase-retrieval dictionary object and descriptions
PH_RETRIEVAL = {
    'h':phretrieve.hilbert_retrieval,
    'w':phretrieve.wave_retrieval,
    'c':phretrieve.correlation_retrieval,
    'w2':phretrieve.wave_retrieval_2,
    'cw':phretrieve.compressed_wave_retrieval,
    'ho':phretrieve.holographic_retrieval,
    'cho':phretrieve.compressed_holographic_retrieval,
    't':phretrieve.test_retrieval}
PH_RETRIEVAL_DESC = {
    'h':'Hilbert-transform phase retrieval method',
    'w':'Gaussian filter over k-region of choice',
    'c':'Correlation method of retrieval',
    't':'Test method of retrieval'}

# operations as python dictionary objects
# one dictionary contains function references, 
# the other contains function descriptions
OPERATIONS = {
    'p':ph_retrieve,
    's':ph_steps,
    'm':ph_movie,
    't':ph_temporal,
    'q':ph_quit}
OPERATION_DESC = {
    'p':'Retrieves phase from selected frame',
    's':'Plots the steps along the phase-processing operation',
    'm':'Makes a movie of the phase images',
    't':'Plotting and analysis of temporal data',
    'q':'Quits the program'}

# a place for all the default settings (for running in non-interactive test mode)
class Defaults(object):
  # operation
  o = 'm'
  # retrieval method
  r = 't'
  # default test folder
  pathext = "res_v_time_Vconst_20000_mV_20151015_152846/"
  # imshow data limits
  clim = [-1.,1.]
  # frame reference region
  refregion = (slice(-30,-20),slice(-30,-20))
  # region to integrate
  integrate = (slice(5,-5),slice(150,152))

if __name__== "__main__":
  INTERACTIVE = False

  print "Welcome to the phase post-processing module designed by Luke Taylor."
  #HOMEPATH = os.environ['HOME']
  HOMEPATH = os.path.expanduser("~")
  print "Below is your home directory, right?\n{0}".format(HOMEPATH)
  if os.path.isdir(HOMEPATH + "/taylo589_2"):
    TAYLO589 = "/taylo589_2"
  elif os.path.isdir(HOMEPATH + "/taylo589orca"):
    TAYLO589 = "/taylo589orca"
  elif os.path.isdir(HOMEPATH + "/taylo589spideroak"):
    TAYLO589 = "/taylo589spideroak"
  else:
    print "ERROR: main path not found. make sure either ~/taylo589_2 or ~/taylo589orca exist and contain desired datasets."
  
  TESTPATH = HOMEPATH + TAYLO589 + "/python/IDTCam/test/testdata/"
  SAMPLEPATH = HOMEPATH + TAYLO589 + "/python/IDTCam/test/testdata/reference_field/reffield.tif"
  if INTERACTIVE:
    # go to the main program loop
    main_loop(TESTPATH,SAMPLEPATH)
  else:
    Dataset = Frameset(TESTPATH + Defaults.pathext,SAMPLEPATH=SAMPLEPATH,\
        retrieval_method=PH_RETRIEVAL[Defaults.r])
    OPERATIONS[Defaults.o](Dataset)
    
  print "EOL"
