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
import ntpath

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

def extract_sources(in_data, ref, RA, DEC, ang_diam):
	"""
	Find sources that are positioned with the dimensions of the image specified.
	
	:param in_data: data of all sources
	:param ref: reference dictionary for column names
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	"""
	sources_data = []
	for ii in range(0,len(in_data)):
		in_ra = float(in_data[ii][ref["RAJ2000"]]); in_dec = float(in_data[ii][ref["DECJ2000"]])
		if (in_ra >= (RA-1.0) and in_ra < (RA+1.0) and in_dec >= (DEC-1.0) and in_dec < (DEC+1.0)):
			sources_data.append(in_data[ii])
	
	return sources_data
	
def to_Aegean_table(source_data, ref, c_freq, base_name, path):
	"""
	configures source data into a format that Aegean/AeRes can understand.
	
	:param source_data: data arrays for sources constrained by RA and DEC values.
	:param ref: reference dictionary for column names
	:param c_freq: central frequency of the data
	"param base_name: The base name for the output Aegean formatted table
	"param path: The path only of the input filename
	"""
	
	out_file = open(path+'\\'+base_name+'_'+c_freq+'_sources.csv','w')
	out_file.write('island,source,background,local_rms,ra_str,dec_str,ra,err_ra,dec,err_dec,peak_flux,err_peak_flux,int_flux,err_int_flux,a,err_a,b,err_b,pa,err_pa,flags,residual_mean,residual_std,uuid,psf_a,psf_b,psf_pa\n')
	
	num_sources = len(source_data)
	for ii in range(0,num_sources):
		out_file.write(str(ii+1)+',') #island: 1 -> num_sources
		out_file.write('0'+',') #source = 0
		for jj in range(ref['background_'+c_freq],ref['local_rms_'+c_freq]+1): # columns encompassing background/rms information
			out_file.write(str(source_data[ii][jj])+',')
		
		for jj in range(ref['ra_str'],ref['err_DECJ2000']+1): # columns encompassing RA and DEC information
			out_file.write(str(source_data[ii][jj])+',')
		
		# Take int_flux <=> peak_flux
		out_file.write(str(source_data[ii][ref['int_flux_'+c_freq]])+','+str(source_data[ii][ref['err_fit_flux_'+c_freq]])+',') #peak_flux + err
		out_file.write(str(source_data[ii][ref['int_flux_'+c_freq]])+','+str(source_data[ii][ref['err_fit_flux_'+c_freq]])+',') #int_flux + err
	
		# write Gaussian parameters
		out_file.write(str(source_data[ii][ref['a_'+c_freq]])+','+'0.00'+',') # a + err
		out_file.write(str(source_data[ii][ref['b_'+c_freq]])+','+'0.00'+',') # b + err
		out_file.write(str(source_data[ii][ref['pa_'+c_freq]])+','+'0.00'+',') # pa + err
		
		# write flags
		out_file.write(str(source_data[ii][ref['flags_deep']])+',')
		
		#write residual info
		for jj in range(ref['residual_mean_'+c_freq],ref['residual_std_'+c_freq]+1): # residual mean and residual std columns
			out_file.write(str(source_data[ii][jj])+',')
		
		# uuid
		out_file.write('None'+',')
		
		# psf a/b and pa
		for jj in range(ref['psf_a_'+c_freq],ref['psf_pa_'+c_freq]+1): # residual mean and residual std columns
			if jj == ref['psf_pa_'+c_freq]:  # last column -> insert newline char
				out_file.write(str(source_data[ii][jj])+'\n')
			else:
				out_file.write(str(source_data[ii][jj])+',')
	
	out_file.close()
	
def run_AeRes(path, base, fits_filename, c_freq):
	"""
	Runs AeRes.py source subtracting program by Paul Hancock. Outputs a .fits image of the original image with the specified sources subtracted
	
	:param path: the path to the input source table
	:param base: the base name of the file
	:param fits_filename: the file name of the fits file to have sources subtracted from
	:param c_freq: the central frequency.
	"""
	
	if fits_filename == None:
		fits_filename = input(">> Enter .fits file name: ")
	
	os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\AeRes.py'+'"' + ' -c ' + '"'+path+'\\'+base+'_'+c_freq+'_sources.csv'+'"' + ' -f ' + '"'+fits_filename+'"' + ' -r ' + '"'+path+'\\'+base+'_'+c_freq+'_residual.fits'+'"')
	
def main():
	usage = "usage: %prog [options] filename.fits"
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--fitsfile", 
					  action="store",type="string",dest="fits_filename",default=None,
					  help=".fits file name", metavar="FITS_FILE")
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
					  action="store", type="float", dest="ra_map",default=None,
					  help="right ascension of the image",metavar="RA")
	parser.add_option("-d","--dec",
					  action="store", type="float", dest="dec_map",default=None,
					  help="declination of the image",metavar="DEC")
	parser.add_option("-a","--angulardiameter",
					  action="store", type="float", dest="ang_diameter",default=2.0,
					  help="angular diameter of the sides of the image",metavar="ANGULAR_DIAMETER")
	parser.add_option("-b","--base",
					  action="store", type="string", dest="base_name",default="ReSId_out",
					  help="The base name for the output Aegean formatted table")
	
	
	#parser.add_option("-","--",
	#				  action="", dest="",default=,
	#				  help="")
	
	
	(options, args) = parser.parse_args()	
	head, tail = ntpath.split(options.data_filename)
	
	
	# read data from input table
	(in_data, column_ref) = read_data(options.data_filename)
	
	# Extract sources which are constrained by input RA/DEC and ang_diam
	source_data = extract_sources(in_data,column_ref,options.ra_map,options.dec_map,options.ang_diameter)
	
	# Convert source data to Aegean format table
	to_Aegean_table(source_data,column_ref,options.central_freq,options.base_name,head)
	
	run_AeRes(head, options.base_name, options.fits_filename, options.central_freq)
	
	
	
main()