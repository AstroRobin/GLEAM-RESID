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

from progressbar import ProgressBar, Bar, Percentage

def find_filename(filename):
	"""
	Find .fits file in directory.
	
	:param filename: The name of the fits file to be searched for
	
	:return: The .fits file name and path
	"""
	if (verbose): print "\n <Searching for .fits file>\n  ** Searching for '"+filename+"' **"
	
	found_filenames = []
	dirs = ['CDFS','ELAIS_S1','COSMOS']
	for ii in range(0,len(dirs)):
		dir_str = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\'+dirs[ii]
		for dir_filename in os.listdir(dir_str):
			if ".fits" in dir_filename:
				if filename in dir_filename:
					found_filenames.append(dirs[ii]+'\\'+dir_filename)
	if (len(found_filenames) == 1):
		print " ** Found .fits file: ", found_filenames[0], " **"
		filename = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\' + found_filenames[0]
	elif (len(found_filenames) > 1):
		print "  ** .fits files found: ** "
		for kk in range(0,len(found_filenames)):
			print " ", kk+1, " - ",found_filenames[kk]
			file_choice = "-1"
		while (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
			file_choice = input(">> Select file: ")
			if (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
				print " ** invalid choice ** "
		filename = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR3\\' + found_filenames[int(file_choice)-1]
	else:
		print " ** No .fits files found with name '", filename,"' **\n   -- ABORTING --   "
		exit()
		
	return filename
	
def find_gal_filename(filename):
	"""
	Find .fits file in directory for a dwarf galaxy
	
	:param filename: The name of the fits file to be searched for
	
	:return: The .fits file name and path
	"""
	if (verbose): print "\n <Searching for .fits file>\n  ** Searching for '"+filename+"' **"

	found_filenames = []
	dir_str = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images'
	gal_dirs = os.listdir(dir_str)
	for dir_filename in gal_dirs:
		if (filename in dir_filename):
			found_filenames.append(dir_filename)
	if (len(found_filenames) == 1):
		print "  ** Found .fits file: ", found_filenames[0], " **"
		filename = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\'+found_filenames[0]+'\\'+found_filenames[0]+'.fits'
	elif (len(found_filenames) > 1):
		print "  ** .fits files found: ** "
		for kk in range(0,len(found_filenames)):
			print " ", kk+1, " - ",found_filenames[kk]
		file_choice = "-1"
		while (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
			file_choice = input(">> Select file: ")
			if (int(file_choice) < 1 or int(file_choice) > len(found_filenames)):
				print " ** invalid choice ** "
		filename = 'C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\'+found_filenames[int(file_choice)-1]+'\\'+found_filenames[int(file_choice)-1]+'.fits'
	else:
		print " ** No .fits files found with name '", filename,"' **\n   -- ABORTING --   "
		exit()
		
	return filename
	

def read_data(filename, RA, DEC, ang_diam, head):
	"""
	Read data from file.
	Take input data of sources and convert to a table
	
	data format: [[col1,col2,col3,...],[col1,col2,col3,...],[col1,col2,col3,...],...,n_rows]
	
	dict format: {0:Names,1:Background_deep,...,34:int_flux_84,...,n_columns}
	
	:param filename: name of the input filename
	:param RA: right ascension of the map
	:param DEC: declination of the map
	:param ang_diam: angular diameter of the fits image
	:param head: path for where to search for pre existing data tables
	
	:return: The data given from the input table in an astropy Table
	"""
	if (verbose): print "\n <Reading in data> \n  ** Searching for pre-existing data file **"
	
	catch = False
	# Searching for filenames of form "GLEAM_chunk_{RA}_{DEC}_{ang_diam}.fits" and encompass input position.
	for search_filename in os.listdir(head):
		if (catch == True):
			break
		if ('.fits' in search_filename and 'GLEAM_chunk' in search_filename):	
			(RA_file,DEC_file,ang_diam_file) = search_filename.replace('.fits','').split('_')[2:5]
			RA_file = float(RA_file); DEC_file = float(DEC_file); ang_diam_file = float(ang_diam_file)
			if (RA - ang_diam >= RA_file - ang_diam_file and RA + ang_diam <= RA_file + ang_diam_file and DEC - ang_diam >= DEC_file - ang_diam_file and DEC + ang_diam <= DEC_file + ang_diam_file):
				print "\n  ** Found pre-existing data file: ", search_filename, " **"
				print "     - RA: ", RA_file, "\n     - DEC: ", DEC_file, "\n     - Angular diameter: ", ang_diam_file
				print "\n       Use this file? "
				choice = '' # enter while loop
				while (choice != 'y' and choice != 'n' and choice != 'Y' and choice != 'N'):
					choice = input(">> (y/n)?: ")
					if (choice == 'y' or choice == 'Y'):
						filename = head + "\\" + search_filename
						catch = True
					elif (choice == 'n' or choice == 'N'):
						print " ** Searching for other files ** "
					else:
						print " ** invalid choice ** "
	if (catch == False):
		print " ** No files pre-existing with appropriate positional parameters ** "
		
	if verbose: print "\n  ** Using input data file: **\n "+ filename
	in_data = Table.read(filename)
	
	return in_data
	

def extract_sources(data, RA, DEC, ang_diam, head):
	"""
	Find sources that are positioned with the dimensions of the image specified.
	
	:param in_data: data of all sources
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	:param head: path to destination folder
	
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
	
	ang_diam_buff = ang_diam*math.sqrt(2)
	num_sources = len(data)
	
	found_sources = 0
	pbar = ProgressBar(widgets=['  ** Extracting sources: ', Percentage(), ' ', Bar(), ' ** '],maxval=num_sources).start()
	for ii in range(0,num_sources):
		pbar.update(ii+1)
		if (data[ii]['RAJ2000'] >= RA - 0.5*ang_diam_buff and data[ii]['RAJ2000'] <= RA + 0.5*ang_diam_buff and data[ii]['DECJ2000'] >= DEC - 0.5*ang_diam_buff and data[ii]['DECJ2000'] <= DEC + 0.5*ang_diam_buff):
			source_data.add_row(data[ii])
			found_sources += 1
	pbar.finish()
	
	print "\n  ** Position bounds: \n      - RA: ",RA - 0.5*ang_diam_buff," -> ",RA + 0.5*ang_diam_buff," \n      - DEC: ",DEC - 0.5*ang_diam_buff," -> ",DEC + 0.5*ang_diam_buff, "\n  ** Number of sources sources found: ", len(source_data)
	
	catch = False
	filename = "GLEAM_chunk_"+str(RA)+"_"+str(DEC)+"_"+str(ang_diam)+".fits"
	for search_filename in os.listdir(head):
		if (catch): break
		if (filename == search_filename):
			print "\n  ** File already exists: **  \n ", filename
			catch = True
	if (catch == False):
		print "\n  ** Writing to file: **  \n ", filename
		source_data.write(head+"\\"+filename)
	
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
	
	
def to_Aegean_table(in_data, c_freq, RA, DEC, ang_diam, head):
	"""
	configures source data into a format that Aegean/AeRes can understand.
	
	:param source_data: data arrays for sources constrained by RA and DEC values.
	:param c_freq: central frequency of the data
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	"param head: The path for the data 'snippet'; i.e. the Aegean formatted single frequency column table of the in_data
	"""
	
	num_sources = len(in_data)
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
	for ii in range(0,num_sources):
		(peak_flux_arr[ii], err_peak_flux_arr[ii]) = calc_peak_flux(in_data['a_'+c_freq][ii],in_data['b_'+c_freq][ii],in_data['psf_a_'+c_freq][ii],in_data['psf_b_'+c_freq][ii],in_data['int_flux_'+c_freq][ii],in_data['err_fit_flux_'+c_freq][ii])
	out_data['peak_flux'] = peak_flux_arr
	out_data['err_peak_flux'] = err_peak_flux_arr
	out_data['int_flux'] = in_data['int_flux_'+c_freq]
	out_data['err_int_flux'] = in_data['err_fit_flux_'+c_freq]
	
	# Source shape information
	out_data['a'] = in_data['a_'+c_freq]
	out_data['err_a'] = 0.0 # no a_err given in GLEAMIDR3.fits
	out_data['b'] = in_data['b_'+c_freq]
	out_data['err_b'] = 0.0 # no b_err given in GLEAMIDR3.fits
	out_data['pa'] = in_data['pa_'+c_freq]
	out_data['err_pa'] = 0.0 # no pa_err given in GLEAMIDR3.fits
	
	# Flags information
	out_data['flags'] = in_data['flags_deep']
	
	# Residuals information
	out_data['residual_mean'] = in_data['residual_mean_'+c_freq]
	out_data['residual_std'] = in_data['residual_std_'+c_freq]
	
	# uuid information
	out_data['uuid'] = ['None']*num_sources
	
	# psf information
	out_data['psf_a'] = in_data['psf_a_'+c_freq]
	out_data['psf_b'] = in_data['psf_b_'+c_freq]
	out_data['psf_pa'] = in_data['psf_pa_'+c_freq]
	
	# edit this output name to something appropriate
	filename = head+"\\"+"GLEAM_snippet_"+str(RA)+"_"+str(DEC)+"_"+str(ang_diam)+'_'+c_freq+'.fits'
	print  filename
	out_data.write(filename,format='fits')

	
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
	parser.add_option("-v","--verbose",
					  action="store_true", dest="verbose",default=False,
					  help="print status messages to stdout")
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
	global verbose; verbose = options.verbose
	
	if (options.galaxy_name != None and options.fits_filename != None):
		print " ** -g (--galaxy) and -f (--fitsfile) have been specified **\n   -- ABORTING --   "
		exit()
	
	if (options.galaxy_name != None):
		fits_filename = find_gal_filename(options.galaxy_name)	
	else:
		if (options.fits_filename != None):
			fits_filename = options.fits_filename
		else:
			fits_filename = input(' >> Search for .fits file: ')
		fits_filename = find_filename(fits_filename)

	head, tail = ntpath.split(fits_filename)
	if (verbose): print "\n  ** Using .fits file name: '.../"+tail+"' ** "
	
	# read data from input table
	in_data = read_data(options.data_filename,float(options.ra_map),float(options.dec_map),float(options.ang_diameter),head)
	
	# Extract sources which are constrained by input RA/DEC and ang_diam
	source_data = extract_sources(in_data,options.ra_map,options.dec_map,options.ang_diameter,head)
	
	# Convert source data to Aegean format table
	if (options.base_name == None):
		base = tail.replace('"','').replace('.fits','')
	else:
		base = options.base_name
	
	to_Aegean_table(source_data,options.central_freq,options.ra_map,options.dec_map,options.ang_diameter,head)
	
	exit()
	run_AeRes(head, base, fits_filename, options.central_freq)
	
main()
