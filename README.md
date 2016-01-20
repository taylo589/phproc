# phproc
Image Processing for Phase Retrieval from Interferograms

This project is intended as a python module for processing interferograms to retrieve optical phase information. 

Image processing is a technique for extracting information from images (there are probably better definitions).

Interferograms are images which display an interference pattern, generally from known sources of light.

Based on a priori knowledge (i.e. from whomever set up the optical experiment used to generate the interferograms), previously-invisible information about objects may be retrieved.

Optical phase information is information about objects and materials that gets recorded in an optical wave when the wave interacts (reflects or transmitts) through the object.

This module serves as a short collection of techniques for retrieving the optical phase from interferometric data frames.

See my paper at [http://dx.doi.org/10.1364/AO.54.009010]

Also see the image below for a description of the Hilbert phase processing algorithm:

![Phase-processing algorithm](https://github.com/taylo589/phproc/blob/master/Figure2.png "Phase-processing Algorithm")
