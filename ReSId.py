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
	
def read_data(filename):
	file = open(filename,'r')
	column_names = file.readline().replace('\n','').split(','); num_columns = len(column_names); #*print column_names
	lines = file.readlines()
	
	in_data = []
	for ii in range(0,len(lines)):
		in_data.append(lines[ii].replace('\n','').split(','))
	
	print in_data
		
def main():
	usage = "usage: %prog [options] filename.fits"
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="write report to FILE", metavar="FILE")
	parser.add_option("-q", "--quiet",
					  action="store_false", dest="verbose", default=True,
					  help="don't print status messages to stdout")
	parser.add_option("-c", "--centralfreq", 
					  action="store",type="string",dest="central_freq",
					  help="provide central frequency", metavar="FREQ")
	parser.add_option("-d","--datafile",
					  action="store", dest="data_filename", 
					  default="C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\GLEAMIDR3.csv",
					  help="destination of input table for sources",metavar="SOURCESFILE")
	
	(options, args) = parser.parse_args()
	
	read_data(options.data_filename)
	

main()