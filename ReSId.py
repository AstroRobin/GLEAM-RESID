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
	"""
	Read data from file.
	Take input data of sources and separate into rows that have each been split by column. Also create a dictionary reference of the column names.
	
	data format: [[col1,col2,col3,...],[col1,col2,col3,...],[col1,col2,col3,...],...,n_rows]
	
	dict format: {0:Names,1:Background_deep,...,34:int_flux_84,...,n_columns}
	
	:param filename: name of the input filename
	:return: The data given from the table in array form, plus a dictionary reference of the column names
	"""

	file = open(filename,'r')
	column_names = file.readline().replace('\n','').split(','); num_columns = len(column_names); #*print column_names
	lines = file.readlines()
	
	in_data = []
	for ii in range(0,len(lines)):
		in_data.append(lines[ii].replace('\n','').split(','))
	
	column_ref = {};
	for ii in range(0,num_columns):
		column_ref[column_names[ii]] = ii
	
	return (in_data,column_ref)

def extract_sources(in_data, RA, DEC, ang_diam):
	"""
	Find sources that are positioned with the dimensions of the image specified.
	
	:param in_data: data of all sources
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	"""
		
	
def main():
	usage = "usage: %prog [options] filename.fits"
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--file", dest="filename",
					  help="write report to FILE", metavar="FILE")
	parser.add_option("-q", "--quiet",
					  action="store_false", dest="verbose", default=True,
					  help="don't print status messages to stdout")
	parser.add_option("-c", "--centralfreq", 
					  action="store",type="string",dest="central_freq",
					  help="provide central frequency", metavar="FREQUENCY")
	parser.add_option("-i","--datafile",
					  action="store", dest="data_filename", 
					  default="C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\CDFS\\CDFS_sources.csv",
					  help="destination of input table for sources",metavar="SOURCES_FILE")
	parser.add_option("-r","--ra",
					  action="store", dest="ra_map",default=None,
					  help="right ascension of the image",metavar="RA")
	parser.add_option("-d","--dec",
					  action="store", dest="dec_map",default=None,
					  help="declination of the image",metavar="DEC")
	parser.add_option("-a","--angulardiameter",
					  action="store", type="double", dest="ang_diameter",default=2.0,
					  help="angular diameter of the sides of the image",metavar="ANGULAR_DIAMETER")
	
	
	#parser.add_option("-","--",
	#				  action="", dest="",default=,
	#				  help="")
	

	(options, args) = parser.parse_args()	
	
	(in_data, column_ref) = read_data(options.data_filename)
	
	extract_sources(in_data,options.RA,options.DEC,options.ang_dimaeter)


main()