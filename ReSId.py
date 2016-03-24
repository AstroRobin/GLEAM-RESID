"""
The GLEAM REsidual Source IDentifier program

Created by:
Robin Cook
March 24 2016

Modifications by:
Robin Cook
"""

# Imports
import sys
import os
import numpy as np
import math

# Other Imports
import scipy
import astropy

from optparse import OptionParser

usage = "usage: %prog [options] filename.fits"
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option("-c", "--freq", 
				  action="store",type="string",dest="central_freq")
				  #,help="provide central frequency", metavar="FREQ")

(options, args) = parser.parse_args()

print options.central_freq
