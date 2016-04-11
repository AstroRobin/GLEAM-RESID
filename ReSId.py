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
from astropy.table import Table
import ntpath

from optparse import OptionParser

def find_filename(filename):
	"""
	Find .fits file in directory.
	
	:param filename: name of the file being searched for
	
	:return: The .fits filename and path
	"""
	
	found_filenames = []
	
	dirs = ['CDFS','ELAIS_S1','COSMOS']
	for ii in range(0,len(dirs)):
		dir_str = "C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\"+dirs[ii]
		for dir_filename in os.listdir(dir_str):
			if ".fits" in dir_filename:
				if filename in dir_filename:
					found_filenames.append(dirs[ii]+'\\'+dir_filename)	
	if (len(found_filenames) == 1):
		print " ** Found .fits file: ", found_filenames[0], " **"
		filename = '"C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\"' + found_filenames[0] + '"'
	elif (len(found_filenames) > 1):
		for kk in range(0,len(found_filenames)):
			print kk+1, " - ",found_filenames[kk]
			file_choice = "-1"
		while (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
			file_choice = input(">> Select file: ")
			print int(file_choice)
			if (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
				print " ** invalid choice ** "
		filename = '"C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\' + found_filenames[int(file_choice)-1] + '"'
	else:
		print " ** No .fits files found with name '", filename,"' **\n   -- ABORTING --   "
		exit()
		
	return filename
	
def find_gal_filename(filename):
	"""
	Find .fits file in directory for a dwarf galaxy
	
	:return: The .fits file name and path
	"""
	
	
	found_filenames = []
	dir_str = "C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images"
	gal_dirs = os.listdir(dir_str)
	for dir_filename in gal_dirs:
		if (filename in dir_filename):
			found_filenames.append(dir_filename)
	if (len(found_filenames) == 1):
		print " ** Found .fits file: ", found_filenames[0], " **"
		filename = '"C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\'+found_filenames[0]+'\\'+found_filenames[0]+'.fits"'
	elif (len(found_filenames) > 1):
		for kk in range(0,len(found_filenames)):
			print kk+1, " - ",found_filenames[kk]
			file_choice = "-1"
		while (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
			file_choice = input(">> Select file: ")
			print int(file_choice)
			if (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
				print " ** invalid choice ** "
		filename = '"C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\'+found_filenames[int(file_choice)-1]+'\\'+found_filenames[int(file_choice)-1]+'.fits"'
	else:
		print " ** No .fits files found with name '", filename,"' **\n   -- ABORTING --   "
		exit()
		
	return filename
	
	
	
def read_data(filename, path):
	"""
	Read data from file.
	Take input data of sources and convert to a table
	
	data format: [[col1,col2,col3,...],[col1,col2,col3,...],[col1,col2,col3,...],...,n_rows]
	
	dict format: {0:Names,1:Background_deep,...,34:int_flux_84,...,n_columns}
	
	:param filename: name of the input filename
	:param path: path for where to search for pre existing tables
	:return: The data given from the input table in an astropy Table
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

def extract_sources(data, RA, DEC, ang_diam):
	"""
	Find sources that are positioned with the dimensions of the image specified.
	
	:param in_data: data of all sources
	:param ref: reference dictionary for column names #edit
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	
	:return: a table of sources which were found to lie within the specified coordinates
	"""
	
	dt = []
	for ii in data.colnames:
		if ('Name' in ii or '_str' in ii):
			dt.append('S20')
		elif('flags' in ii):
			dt.append('int')
		else:
			dt.append('float')
			
	source_data = Table(names=(data.colnames),dtype=dt)
	
	ang_diam = ang_diam*math.sqrt(2)
	
	for ii in range(0,len(data_table)):
		if (data[ii]['RAJ2000'] >= RA - 0.5*ang_diam and data[ii]['RAJ2000'] <= RA + 0.5*ang_diam and data[ii]['DECJ2000'] >= DEC - 0.5*ang_diam and data[ii]['DECJ2000'] <= DEC + 0.5*ang_diam):
			source_data.add_row(data[ii])
	
	#sources_data = []
	#for ii in range(0,len(in_data)):
	#	in_ra = float(in_data[ii][ref["RAJ2000"]]); in_dec = float(in_data[ii][ref["DECJ2000"]])
	#	if (in_ra >= (RA-0.5*ang_diam) and in_ra < (RA+0.5*ang_diam) and in_dec >= (DEC-0.5*ang_diam) and in_dec < (DEC+0.5*ang_diam)):
	#		sources_data.append(in_data[ii])
	#return sources_data
	
	print " ** Position bounds: \n   - RA: ",RA - 0.5*ang_diam," -> ",RA + 0.5*ang_diam," \n   - DEC: ",DEC - 0.5*ang_diam," -> ",DEC - 0.5*ang_diam, "\n  ** Number of sources sources found: ", len(source_data)
	
	return source_data
	
def calc_peak_flux (a,b,psf_a,psf_b,int_flux,err_int_flux):
	"""
	Calculates th peak flux of the source using the equation: peak_flux = int_flux x (psf_a x psf_b)/(a x b) 
	
	:param a: source FWHM semi-major axis
	:param b: source FWHM semi-minor axis
	:param psf_a: synthesized beam FWHM semi-major axis
	:param psf_a: synthesized beam FWHM semi-major axis
	:param int_flux: integrated flux of the source
	:param err_int_flux: error in the integrated flux of the source
	
	:return: peak_flux and err_peak_flux
	"""
	
	peak_flux = ((psf_a*psf_b)/(a*b))*int_flux
	err_peak_flux = (peak_flux/int_flux)*err_int_flux
	
	return [peak_flux, err_peak_flux]
	
	
def to_Aegean_table(in_data, ref, c_freq, base_name,path):
	"""
	configures source data into a format that Aegean/AeRes can understand.
	
	:param source_data: data arrays for sources constrained by RA and DEC values.
	:param ref: reference dictionary for column names
	:param c_freq: central frequency of the data
	"param base_name: The base name for the output Aegean formatted table
	"""
	
	#out_file = open(path+'\\'+base_name+'_'+c_freq+'_sources.csv','w')
	#out_file.write('island,source,background,local_rms,ra_str,dec_str,ra,err_ra,dec,err_dec,peak_flux,err_peak_flux,int_flux,err_int_flux,a,err_a,b,err_b,pa,err_pa,flags,residual_mean,residual_std,uuid,psf_a,psf_b,psf_pa\n')
	
	num_sources = len(source_data)
	out_data = Table()
	
	# Aegean requirement information
	out_data['island'] = [kk for kk in range(1,num_sources+1)] #island: 1 -> num_sources
	out_data['source'] = [0]*num_sources # source -> 0
	
	# Background rms information
	out_data['background'] = in_data['background_'+c_freq] # background noise for this frequency band
	out_data['local_rms'] = in_data['local_rms_'+c_freq] # local rms for this frequency band
	
	# Positional information
	out_data['ra_str'] = in_data['ra_str']; out_data['dec_str'] = in_data['dec_str'] # RA/DEC Strings
	out_data['ra'] = in_data['RAJ2000']; out_data['err_ra'] = in_data['err_RAJ2000'] # RA + RA_err
	out_data['dec'] = in_data['DECJ2000']; out_data['err_dec'] = in_data['err_DECJ2000'] # DEC + DEC_err
	
	# Peak and integrated flux data
	peak_flux_arr = [0]*num_sources; err_peak_flux_arr = [0]*num_sources
	for ii in range(0,num_source):
		(peak_flux_arr[ii], err_peak_flux_arr[ii]) = calc_peak_flux(in_data['a_'+c_freq][ii],in_data['b_'+c_freq][ii],in_data['psf_a_'+c_freq][ii],in_data['psf_b_'+c_freq][ii],in_data['int_flux_'+c_freq][ii],in_data['err_fit_flux_'+c_freq][ii])
	out_data['peak_flux'] = peak_flux_arr
	out_data['err_peak_flux'] = err_peak_flux_arr
	out_data['int_flux'] = in_data['int_flux_'+c_freq]
	out_data['err_int_flux'] = in_data['err_fit_flux_'+c_freq]
	
	# Source shape information
	for ii in ['a','err_a','b','err_b','pa','err_pa']:
		out_data[ii] = in_data[ii+'_'+c_freq]
	
	# Flags information
	out_data['flags'] = in_data['flags_deep']
	
	# Residuals information
	out_data['residual_mean'] = in_data['residual_mean_'+c_freq]
	out_data['residual_std'] = in_data['residual_std_'+c_freq]
	
	# uuid information
	out_data['uuid'] = [None]*num_sources
	
	# psf information
	out_data['psf_a'] = in_data['psf_a_'+c_freq]
	out_data['psf_b'] = in_data['psf_b_'+c_freq]
	out_data['psf_pa'] = in_data['psf_pa_'+c_freq]
	
	out_data.write(head+'\\'+base+'_resid.fits')
	
def run_AeRes(path, base, fits_filename, c_freq):
	"""
	Runs AeRes.py source subtracting program by Paul Hancock. Outputs a .fits image of the original image with the specified sources subtracted
	
	:param path: the path to the input source table
	:param base: the base name of the file
	:param fits_filename: the file name of the fits file to have sources subtracted from
	:param c_freq: the central frequency.
	"""
	
	os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\AeRes.py'+'"' + ' -c ' + '"'+path+'\\'+base+'_'+c_freq+'_sources.csv'+'"' + ' -f ' + '"'+fits_filename+'"' + ' -r ' + '"'+path+'\\'+base+'_'+c_freq+'_residual.fits'+'"')

def main():
	usage = "usage: %prog [options] filename.fits"
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--fitsfile", 
					  action="store",type="string",dest="fits_filename",default=None,
					  help=".fits file name", metavar="FITS_FILE")
	parser.add_option("-g","--galaxy",
					  action="store", type="string", dest="galaxy_name",default=None,
					  help="The name of the Dwarf galaxy",metavar="GALAXY_NAME")
	parser.add_option("-q", "--quiet",
					  action="store_false", dest="verbose", default=True,
					  help="don't print status messages to stdout")
	parser.add_option("-c", "--centralfreq", 
					  action="store",type="string",dest="central_freq",
					  help="provide central frequency", metavar="FREQUENCY")
	parser.add_option("-i","--datafile",
					  action="store", dest="data_filename", 
					  default="C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\GLEAMIDR3.csv",
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
					  action="store", type="string", dest="base_name",default=None,
					  help="The base name for the output Aegean formatted table")
	
	
	#parser.add_option("-","--",
	#				  action="", dest="",default=,
	#				  help="")
	
	
	(options, args) = parser.parse_args()	
	
	if (options.galaxy_name != None):
		# this will not work
		fits_filename = find_gal_filename(options.galaxy_name)	
	else:
		if (options.fits_filename != None):
			fits_filename = options.fits_filename
		else:
			fits_filename = input(' >> Search for .fits file: ')
		fits_filename = find_filename(fits_filename)
	
	print fits_filename
	
	head, tail = ntpath.split(fits_filename)
	
	
	# read data from input table
	in_data = Table.read(options.data_filename)
	exit()
	# Extract sources which are constrained by input RA/DEC and ang_diam
	source_data = extract_sources(in_data,column_ref,options.ra_map,options.dec_map,options.ang_diameter)
	
	# Convert source data to Aegean format table
	if (options.base_name == None):
		base = tail.replace('"','').replace('.fits','')
	else:
		base = options.base_name
		
	to_Aegean_table(source_data,options.central_freq,base,head)
	
	run_AeRes(head, base, fits_filename, options.central_freq)
	
main()
