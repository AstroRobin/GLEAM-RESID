"""
The GLEAM REsidual Source IDentifier program (ReSId)

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
import time

# Other Imports
import scipy
import astropy
from astropy.table import Table
from astropy import units as u
from astropy.wcs import wcs
import astropy.io.fits as pyfits
import ntpath

from optparse import OptionParser

from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory


# Imports for get_gleam()
# need to get Nick to install wsgi_intercept
#from gleam_vo_example import GleamVoProxy, download_file
#import pyvo

# Hardcoding
IDR_version = "5"

def choose(files):
	"""
	Given a list of filenames, this function will list all filenames in a formatted order and prompt the user to select a file
	
	<param: filenames> - a list of filenames
	
	<return: filename> - the chosen filename
	"""
	
	print "  ** Multiple ({0}) files found  ** ".format(len(files))
	for kk in range(0,len(files)): print " {0} - {1}".format(kk+1,files[kk])
	while True:
			choice = '1' if (auto) else raw_input("\n >> Choose file: ") # auto_answer -> first file in list
			try:
				choice = int(choice)
				if (choice >= 1 and choice <= len(files)):
					return files[choice-1]; break
				else:
					print "  ** ERROR: input out of bounds **  "
			except ValueError:
				print "  ** ERROR: invalid input **  "

def find_file(search_path):
	"""
	Find .fits file in directory.
	
	<param: search_path> - The name of the fits file (and path) to be searched for
	
	<return: file path>
	"""
	
	dir, in_filename = ntpath.split(search_path)
	if (verbose): print "\n <Searching for .fits file>\n  ** Searching for '{0}' in .../{1[0]}/{1[1]} **".format(filename,dir.split("\\")[-2:])

	found_files = []
	for filename in os.listdir(dir):
		if (in_filename in filename):
			found_files.append(filename)
	if (len(found_files) == 1): # if only one appropriate file found
		if (verbose): print "  ** Found .fits file: {0} **".format(found_files[0])
		file = "{0}\\{1}".format(dir,found_files[0])
	elif (len(found_files) > 1): # if multiple appropriate files found
		file = "{0}\\{1}".format(dir,choose(found_files))
	else:  # if no appropriate files found
		if (verbose): print " ** \"{0}\" not found in .../{1[0]}/{1[1]} **".format(in_filename,dir.split("\\")[-2:])
		file = ""
	return file
	
def get_position(ra,dec,ang_diam):
	"""
	get position parameters from user.
	
	<param: ra> - right ascension
	<param: dec> - declination
	<param: ang_diam> -	angular diameter
	
	<return: [ra,dec,ang_diam]> - right ascension, declination and angular diameter after validity checks
	"""
	if (verbose): print "\n <Checking positional parameters>\n"
	
	# Check if RA, DEC, and Angular diameter are given, if not: enter now
	if (ra == None): # no RA given
		while (True):
			try: 
				ra = float(raw_input("\n>> Enter Right Ascension: "))
				if (ra >= 0.0 and ra < 360.0): break
				else: print "\n ** WARNING: RA out of bounds (0 < RA < 360)**"
			except ValueError: print "\n  ** ERROR: invalid input **"
	if (dec == None): # no DEC given 
		while (True):
			try: 
				dec = float(raw_input("\n>> Enter Declination: "))
				if (dec >= -90.0 and dec <= +90.0): break
				else: print "\n ** WARNING: DEC out of bounds (-90 < DEC < +90) **"
			except ValueError: print "\n  ** ERROR: invalid input **"
	if (ang_diam == None): # no angular diameter given
		while (True):
			try: 
				ang_diam = float(raw_input("\n>> Enter Angular Diameter: "))
				if (ang_diam > 0.0 and ang_diam < 5.0): break
				else: print "\n ** WARNING: Angular diameter out of bounds (ang. diameter < 5 degrees) **"
			except ValueError: print "\n  ** ERROR: invalid input **"

	return (ra,dec,ang_diam)		
			
def find_gal_file(gals_dir,galaxy,ra,dec,ang_diam,freq):
	"""
	Firstly, find whether or not the given galaxy name exists within the galaxies directory - if it does not exists, ask the user whether they wish to create a new galaxy directory and set this as the directory to look for files in.
	If the galaxy's directory does exist, set this as the directory to look for files within. Secondly search for a .fits file with the key words 'cutout' and the specified frequency. If it does exist, return the path.
	If the GLEAM cutout file does not exist (as is the case when a new galaxy directory has been created), propmt the user whether they wish to dowload the cutout. If yes, use the RA, DEC, angular diameter and frequency specified in ReSId.py options, else ask the user for RA, DEC and angular diameter; Note, default frequency is 'wide'.
	
	<param: gals_dir> - The directory in which to look for galaxy folders
	<param: galaxy> - The name of the fits file to be searched for
	<param: RA> - right ascension of the map
	<param: DEC> - declination of the map
	<param: ang_diam> - angular diameter of the fits image
	<param: freq> - the frequency of the image
	
	:return: The .fits file name and path
	"""
	if (verbose): print "\n <Searching for .fits file>\n\n  ** Searching for galaxy: '{0}' **".format(galaxy)
	
	# dir = "C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\" # root directory for DSph galaxy images
	# look for directories with the name of the galaxy given
	catch = False
	print "\n'{0}'\n".format(gals_dir)
	for dir_name in os.listdir(gals_dir): # iterate over all galaxy directories
		if (catch): break
		if (dir_name == galaxy):
			if (verbose): print "  ** Found directory: '{0}' **  \n".format(galaxy)
			dir = "{0}/{1}".format(gals_dir,galaxy)
			catch = True; break
	if (catch == False):
		print "\n  ** WARNING: no directory \"{0}\" found in \"{1}\"**\n  ".format(galaxy,gals_dir)
		while True:
			choice = 'y' if (auto) else str(raw_input("\n >> Make new galaxy directory \"{0}\" (y/n)?: ".format(galaxy))) #auto_answer -> yes, create new galaxy directory
			if ('y' in choice.lower()): # make new directory for this folder
				os.system("mkdir \"{0}\"".format(dir+"/"+galaxy))
				if (verbose): print "  ** Directory: '{0}' has been created. **".format(galaxy)
				dir = "{0}/{1}".format(gals_dir,galaxy)
				catch = True; break
			elif ('n' in choice.lower()): # don't make new directory -> ABORT
				if (verbose): print "\n    -- ABORTING --   "; exit()
	
	
	# look for files in galaxy directory with 'cutout' in their name
	if (verbose): print "  ** Searching for appropriate 'GLEAM_cutout' in '.../{0}' ** ".format(dir.split("/")[-1])
	
	found_files = []
	for filename in os.listdir(dir):
		if ("cutout" in filename and freq in filename):
			found_files.append(filename)
	if (len(found_files) == 1): # if only one appropriate file found
		if (verbose): print "  ** Found .fits file: {0} **".format(found_files[0])
		file = "{0}/{1}".format(dir,found_files[0])
	elif (len(found_files) > 1): # if multiple appropriate files found
		print "  ** Multiple ({0}) files found  ** ".format(len(found_files))
		file = "{0}/{1}".format(dir,choose(found_files))
	else: # if no appropriate files found
		print "  ** WARNING: No GLEAM_cutout_.fits files found in '{0}' directory **".format(galaxy)
		while True:
			choice = 'y' if (auto) else str(raw_input(">> Download cutout(y/n)?: ")) # auto_answer -> yes, download cutout
			if ('y' in choice.lower()): # download cutout
				ra, dec, ang_diam = get_position(ra,dec,ang_diam)
				if (verbose): print " ** Downloading GLEAM cutout for \"{0}\" using parameters: **\n    - RA: {1}\n    - DEC: {2}\n    - Angular diameter: {3}\n    - Frequency: {4} MHz".format(galaxy,ra,dec,ang_diam,freq)
				DL_file = get_cutout(ra, dec, freq, ang_diam, download_dir=dir, listf=False)
				file = "{0}/GLEAM_cutout_{1}_{2}.fits".format(dir,freq,galaxy)
				os.system("rename \"{0}\" \"GLEAM_cutout_{1}_{2}.fits\"".format(DL_file,freq,galaxy))
				break
			elif ('n' in choice.lower()): # don't download cutout, return None
				file = ""; break
			else: 
				print "\n  ** ERROR: invalid input **  "
				
	return file
	
def check_for_file(dir, RA, DEC, ang_diam, in_freq="N/A"):
	"""
	Check whether a file already exists in current directory
	
	:param dir: path for where to search for pre existing data tables
	:param RA: right ascension of the map
	:param DEC: declination of the map
	:param ang_diam: angular diameter of the fits image
	:param in_freq: required frequency of the map
	
	:return filename: the filename and path of the found file
	"""

	if (verbose): print "\n <Searching for pre-existing {0} data file> ".format("chunk" if in_freq=="N/A" else "snippet")
	
	# Searching for filenames of form "GLEAM_[chunk/snippet]_{RA}_{DEC}_{ang_diam}_{freq?}.fits" .
	found_filenames = []
	for filename in os.listdir(dir):
		if (in_freq=="N/A"): # user is looking for a GLEAM_chunk
			if (".fits" in filename and "chunk" in filename):
				(RA_file,DEC_file,ang_diam_file) = filename.replace(".fits","").split("_")[2:]
				RA_file = float(RA_file); DEC_file = float(DEC_file); ang_diam_file = float(ang_diam_file)	
				if (RA - ang_diam >= RA_file - ang_diam_file and RA + ang_diam <= RA_file + ang_diam_file and DEC - ang_diam >= DEC_file - ang_diam_file and DEC + ang_diam <= DEC_file + ang_diam_file):
					found_filenames.append(filename)						
		else: # user is looking for a GLEAM_snippet
			if (".fits" in filename and "snippet" in filename):
				(RA_file,DEC_file,ang_diam_file,freq_file) = filename.replace(".fits","").split("_")[2:6]
				RA_file = float(RA_file); DEC_file = float(DEC_file); ang_diam_file = float(ang_diam_file); freq_file = str(freq_file)	
				if (RA - ang_diam >= RA_file - ang_diam_file and RA + ang_diam <= RA_file + ang_diam_file and DEC - ang_diam >= DEC_file - ang_diam_file and DEC + ang_diam <= DEC_file + ang_diam_file and freq_file == in_freq):
					found_filenames.append(filename)
	
	if (len(found_filenames) == 1):
		if (verbose): print "  ** Found pre-existing file '{0}' ** ".format(found_filenames[0])
		return dir + "\\" + found_filenames[0]
	elif (len(found_filenames) > 1):
		print "  ** Multiple ({0}) pre-existing files found  ** ".format(len(found_filenames))
		filename = dir + "\\" + choose(found_filenames)
	else:
		print "  ** WARNING: no appropriate files found - Returning 'None' ** "
		return None
	# run chunk creation process
	
def get_frequency(freq):
	"""
	Standardizes the input frequency to one that ReSId can understand. Accepts many different input formats for frequency.
	
	:param freq: frequency input by the user
	
	:return: 
	freq: for numerical comparison
	"""
	
	if (verbose): print "\n <Finding frequency range>"
	if ('mhz' in freq.lower()): freq = freq.lower().replace('mhz','')
	if (len(freq) == 2): # check if this is a valid 2 digit input, if so add a zero prefix
		try:
			if (int(freq) < 72): print "  ** WARNING: input frequency out of range **\n   -- ABORTING --   "; exit()
			else: freq = '0'+freq
		except ValueError:
			print "  ** WARNING: input frequency not a number **"
	if (freq.lower() in ['red','r']): freq = 'red'
	if (freq.lower() in ['green','g']): freq = 'green'
	if (freq.lower() in ['blue','b']): freq = 'blue'
	if (freq.lower() in ['white','deep','wide','w']): freq = 'wide'
	

	return freq

def log(last_fits_file,last_gal_dir,last_gal_name,last_catalogue_file,last_Aegean_dir):
	"""
	Log usage history information to a file for later reference
	
	<param: last_fits_file> - most recently used .fits file.
	<param: last_gal_dir> - most recently used directory for accessing galaxies.
	<param: last_gal_name> - most recently used galaxy.
	<param: last_catalogue_file> - most recently used directory for accesing the catalogue of sources.
	<param: last_Aegean_dir> - most recently used directory for reference Aegean programs.
	
	<return: None>
	"""
	if (verbose): print "\n <Logging usage to \"ReSId_log.txt\">"
	
	log_file = open("ReSId_log.txt",'w')
	
	log_file.write("########################################\n## Residual Source Identifier (ReSId) ##\n########################################\n")

	# date-time
	log_file.write("\n** last used ReSId.py **\n{0}\n".format(time.strftime("%c")))
	
	# most recents directories
	log_file.write("\n** last .fits file path **\n{0}\n".format(last_fits_file))
	log_file.write("\n** last dSph galaxies directory **\n{0}\n".format(last_gal_dir))
	log_file.write("\n** last dSph galaxy name **\n{0}\n".format(last_gal_name))
	log_file.write("\n** last catalogue file path **\n{0}\n".format(last_catalogue_file))
	log_file.write("\n** last AEGEAN directory **\n{0}\n".format(last_Aegean_dir))
	
	# print #EOF
	log_file.write("\n#EOF")
	
	log_file.close()
	
def read_log():
	"""
	Reads usage history information from ReSId_log.txt
	
	<return: [last_fit,last_gal_dir,last_gal_name,last_catalogue]>
	"""
	
	try:
		file = open("ReSId_log.txt",'r')
		line = file.readline()
		while ("#EOF" not in line):
			line = file.readline()
			if ("last .fits file path" in line): last_fits_file = file.readline().replace('\n','').replace('\r','')
			elif ("last dSph galaxies directory" in line): last_gal_dir = file.readline().replace('\n','').replace('\r','')
			elif ("last dSph galaxy name" in line): last_gal_name = file.readline().replace('\n','').replace('\r','')
			elif ("last catalogue file path" in line): last_catalogue_file = file.readline().replace('\n','').replace('\r','')
			elif ("last AEGEAN directory" in line): last_Aegean_dir = file.readline().replace('\n','').replace('\r','')
			
		return [last_fits_file,last_gal_dir,last_gal_name,last_catalogue_file,last_Aegean_dir]
		
		file.close()
	except IOError:
		if (verbose): print "\n ** No previous ReSId_log.txt file found **"
		return [""]*5
	
def read_data(filename, RA, DEC, ang_diam, head): # ** Now REDUNDANT! **
	
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
	# function should check for existing filename outside of this function, i.e. in main()
	
	if (verbose): print '\n <Reading in data>'
	
	# this needs some cleaning up
	found_filename = check_for_file(head,RA,DEC,ang_diam)
	if (found_filename != None): filename = found_filename
	if verbose: print "\n  ** Using input data file: **\n "+ filename
	in_data = Table.read(filename)
	
	return in_data
	
def extract_sources(catalogue_file, in_RA, in_DEC, Bmaj, Bmin, Bpa, ang_diam, freq, dir, base):
	"""
	Find sources that are positioned with the dimensions of the image specified.
	
	<param: catalogue_file> - path to the catalogue file
	<param: in_RA> - right ascension of the image
	<param: in_DEC> - declination of the image
	<param: Bmaj> - beam major axis
	<param: Bmin> - beam minor
	<param: Bpa> - beam position angle
	<param: ang_diam> - angular diameter of the image
	<param: freq> - the central frequency of the image
	<param: dir> - path to destination folder
	<param: base> - base name for output file
	
	<return> - the name of the catalogue file produced.
	"""
	if (verbose): print "\n <Extracting sources> \n"
	
	# calculate the position bounds
	try:
		RA_min = in_RA - (0.5*ang_diam*math.sqrt(2.0))/(abs(math.cos(math.radians(in_DEC))))
		RA_max = in_RA + (0.5*ang_diam*math.sqrt(2.0))/(abs(math.cos(math.radians(in_DEC))))
	except ZeroDivisionError:
		RA_min = 0.0
		RA_max = 360.0
	DEC_min = in_DEC - 0.5*ang_diam*math.sqrt(2)
	DEC_max = in_DEC + 0.5*ang_diam*math.sqrt(2)
	#except TypeError:
	#	print "\n ** ERROR: No RA or DEC have been specified **\n     - ABORTING -"; exit()
	
	RA_underflow = False; RA_overflow = False
	if (RA_min < 0.0): # i.e. RA overlaps 0h
		RA_min_over = 360.0 + RA_min
		RA_max_over = 360.0
		RA_min = 0.0
		RA_underflow = True		
	elif (RA_max > 360.0): # i.e. RA overlaps 24h
		RA_min_over = 0.0
		RA_max_over = RA_max - 360.0
		RA_max = 360.0
		RA_overflow = True
	else:
		RA_min_over = 0.0
		RA_max_over = 0.0
	
	#stilts_path = "C:\\Users\\user\\Documents\\TopCat\\stilts.jar"
	stilts_path = "/Users/rcook/bin/stilts.jar"
	input_fmt = "fits"
	output_fmt = "csv"
	
	out_file="{0}/GLEAM_catalogue_{1}.{2}".format(dir,base,output_fmt)
	# use stilts tpipe to extract all sources from GLEAM catalogue that fall within RA and DEC constraints + rename columns for Aegean
	if (verbose): print "\n <Using 'stilts tpipe' to extract sources from GLEAMIDR{0}.fits>".format(IDR_version)
	print "\njava -jar {1} tpipe ifmt={2} omode=out out=sources_temp.txt ofmt={4} ...  \"{5}\"\n".format(freq, stilts_path, input_fmt, out_file, output_fmt,catalogue_file)
	os.system("java -jar {1} tpipe ifmt={2} omode=out out=sources_temp.txt ofmt={4} "
			  "cmd='select \"(DEJ2000 > {5} && DEJ2000 < {6} && ((RAJ2000 >= {7} && RAJ2000 <= {8}) || (RAJ2000 >= {9} && RAJ2000 <= {10})))\"' "
			  "cmd='replacecol -name background background_{0} (background_{0})' "
			  "cmd='replacecol -name local_rms local_rms_{0} (local_rms_{0})' "
			  "cmd='replacecol -name ra_str ra_str (ra_str)' "
			  "cmd='replacecol -name dec_str dec_str (dec_str)' "
			  "cmd='replacecol -name ra RAJ2000 (RAJ2000)' "
			  "cmd='replacecol -name err_ra err_RAJ2000 (err_RAJ2000)' "
			  "cmd='replacecol -name dec DEJ2000 (DEJ2000)' "
			  "cmd='replacecol -name err_dec err_DEJ2000 (err_DEJ2000)' "
			  "cmd='replacecol -name peak_flux peak_flux_{0} (peak_flux_{0})' "
			  "cmd='replacecol -name err_peak_flux err_peak_flux_{0} (err_peak_flux_{0})' "
			  "cmd='replacecol -name int_flux int_flux_{0} (int_flux_{0})' "
			  "cmd='replacecol -name err_int_flux err_int_flux_{0} (err_int_flux_{0})' "
			  
			  "cmd='replacecol -name a a_{0} ((({11}/psf_a_{0})*a_{0}))' " # in arcsec
			  #"cmd='replacecol -name a a_{0} (a_{0})' "
			  #"cmd='replacecol -name err_a _{0} (_{0})' "
			  
			  "cmd='replacecol -name b b_{0} ((({12}/psf_b_{0})*b_{0}))' " # in arcsec
			  #"cmd='replacecol -name b b_{0} (b_{0})' "
			  #"cmd='replacecol -name  err_b_{0} (_{0})' "
			  
			  "cmd='replacecol -name pa pa_{0} (pa_{0})' "
			  #"cmd='replacecol -name pa pa_{0} (pa_{0})' "
			  #"cmd='replacecol -name err_pa _{0} (_{0})' "
			  
			  #"cmd='replacecol -name flags _{0} (_{0})' "
			  "cmd='replacecol -name residual_mean residual_mean_{0} (residual_mean_{0})' "
			  "cmd='replacecol -name residual_std residual_std_{0} (residual_std_{0})' "
			  #"cmd='replacecol -name uuid _{0} (_{0})' "
			  "cmd='replacecol -name psf_a psf_a_{0} (psf_a_{0})' "
			  "cmd='replacecol -name psf_b psf_b_{0} (psf_b_{0})' "
			  "cmd='replacecol -name psf_pa psf_pa_{0} (psf_pa_{0})' "
			  "\"{13}\"".format(freq, stilts_path, input_fmt, out_file, output_fmt, DEC_min, DEC_max, RA_min, RA_max, RA_min_over, RA_max_over, Bmaj, Bmin, catalogue_file))
	
	# create a table using Astropy which contains the missing columns in GLEAMIDRn.fits, i.e. island, source, err_a, err_b, err_pa, flags, uuid
	if (verbose): print " <Creating temporary island table>"
	#num_rows = len(Table.read("sources_temp.txt",format='csv'))
	t = Table.read("sources_temp.txt",format='csv')
	num_rows = len(t)
	Table([[kk for kk in range(1,num_rows+1)]]+[[0 for kk in range(num_rows)]]*6,names=("island","source","err_a","err_b","err_pa","flags","uuid")).write("islands_temp.txt",format='csv')
	# +[[(BMaj/t['psf_a_{0}'.format(freq)][kk])*t['a_{0}'.format(freq)][kk] for kk in range(0,num_rows)]]+[[(BMin/t['psf_b_{0}'.format(freq)][kk])*t['b_{0}'.format(freq)][kk] for kk in range(0,num_rows)]]+[[t['pa_{0}'.format(freq)][kk] for kk in range(0,num_rows)]]
	
	# join the GLEAM data and the false columns
	if (verbose): print " <Joining islands table to source data table>"
	os.system("java -jar {0} tjoin nin=2 ifmt1={1} in1=sources_temp.txt ifmt2=csv in2=islands_temp.txt ofmt=csv out=joined_temp.txt".format(stilts_path,output_fmt))
	# output catalogue with Aegean appropriate column names and ordering
	if (verbose): print " <Rearranging columns>"
	os.system("java -jar {0} tpipe ifmt={1} omode=out ofmt={1} out={2} "
			  "cmd='keepcols \"island source background local_rms ra_str dec_str ra err_ra dec err_dec peak_flux err_peak_flux int_flux err_int_flux a err_a b err_b pa err_pa flags residual_mean residual_std uuid psf_a psf_b psf_pa\"' "
			  "joined_temp.txt".format(stilts_path,output_fmt,out_file))
	if (verbose): print " <Deleting temporary files>"
	os.system("rm sources_temp.txt islands_temp.txt joined_temp.txt") # remove temporary files	
		
	if (verbose):
		if (RA_underflow): print "\n  ** Position bounds: \n      - RA: {0} -> {1}\n      - DEC: {2} -> {3}\n  ** Number of sources sources found: {4}".format(RA_min_over,RA_max,DEC_min,DEC_max,num_rows)
		elif (RA_overflow): print "\n  ** Position bounds: \n      - RA: {0} -> {1}\n      - DEC: {2} -> {3}\n  ** Number of sources sources found: {4}".format(RA_min,RA_max_over,DEC_min,DEC_max,num_rows)
		else: print "\n  ** Position bounds: \n      - RA: {0} -> {1}\n      - DEC: {2} -> {3}\n  ** Number of sources sources found: {4}".format(RA_min,RA_max,DEC_min,DEC_max,num_rows)
		
	return out_file
	
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
	# now redundant -> IDR4 has peak fluxes given.
	
	peak_flux = ((psf_a*psf_b)/(a*b))*int_flux
	err_peak_flux = (peak_flux/int_flux)*err_int_flux
	
	return [peak_flux, err_peak_flux]
	
def to_Aegean_table(in_data, c_freq, RA, DEC, ang_diam, head):
	"""
	configures source data into a format that Aegean/AeRes can understand.
	
	<param: source_data> - data arrays for sources constrained by RA and DEC values.
	<param: c_freq> - central frequency of the data
	<param: RA> - right ascension of the image
	<param: DEC> - declination of the image
	<param: ang_diam> - angular diameter of the image
	<param: head> - The path for the data 'snippet'; i.e. the Aegean formatted single frequency column table of the in_data
	
	<return: catalogue_csv_file> - the path to the csv catalogue file
	"""
	# now redundant -> STILTS will rename and format names
	
	num_sources = len(in_data)
	out_data = Table()
	
	# Aegean requirement information
	out_data['island'] = [kk for kk in range(1,num_sources+1)] #island: 1 -> num_sources
	out_data['source'] = [0]*num_sources # source = 0
	
	# Background rms information
	out_data['background'] = in_data['background_'+c_freq] # background noise for this frequency band
	out_data['local_rms'] = in_data['local_rms_'+c_freq] # local rms for this frequency band
	
	# Positional information
	out_data['ra_str'] = in_data['ra_str']; out_data['dec_str'] = in_data['dec_str'] # RA/DEC Strings
	out_data['ra'] = in_data['RAJ2000']; out_data['err_ra'] = in_data['err_RAJ2000'] # RA + RA_err
	out_data['dec'] = in_data['DEJ2000']; out_data['err_dec'] = in_data['err_DEJ2000'] # DEC + DEC_err
	
	# Peak and integrated flux data
	# peak_flux_arr = [0]*num_sources; err_peak_flux_arr = [0]*num_sources
	# for ii in range(0,num_sources):
	#	(peak_flux_arr[ii], err_peak_flux_arr[ii]) = calc_peak_flux(in_data['a_'+c_freq][ii],in_data['b_'+c_freq][ii],in_data['psf_a_'+c_freq][ii],in_data['psf_b_'+c_freq][ii],in_data['int_flux_'+c_freq][ii],in_data['err_fit_flux_'+c_freq][ii])
	out_data['peak_flux'] = in_data['peak_flux_'+c_freq]
	out_data['err_peak_flux'] = in_data['err_peak_flux_'+c_freq]
	out_data['int_flux'] = in_data['int_flux_'+c_freq]
	out_data['err_int_flux'] = in_data['err_int_flux_'+c_freq]
	
	# Source shape information
	out_data['a'] = in_data['a_'+c_freq]
	out_data['err_a'] = 0.0 # no a_err (8GHz) given in GLEAMIDR4.fits
	out_data['b'] = in_data['b_'+c_freq]
	out_data['err_b'] = 0.0 # no b_err (8GHz) given in GLEAMIDR4.fits
	out_data['pa'] = in_data['pa_'+c_freq]
	out_data['err_pa'] = 0.0 # no pa_err (8GHz) given in GLEAMIDR4.fits
	
	# Flags information -> not in IDR4
	out_data['flags'] = 0
	
	# Residuals information
	out_data['residual_mean'] = in_data['residual_mean_'+c_freq]
	out_data['residual_std'] = in_data['residual_std_'+c_freq]
	
	# uuid information
	out_data['uuid'] = ['None']*num_sources
	
	# psf information
	out_data['psf_a'] = in_data['psf_a_'+c_freq]
	out_data['psf_b'] = in_data['psf_b_'+c_freq]
	out_data['psf_pa'] = in_data['psf_pa_'+c_freq]
	
	
	filename = head+"\\"+"GLEAM_snippet_"+str(RA)+"_"+str(DEC)+"_"+str(ang_diam)+'_'+c_freq+'.fits'
	if (verbose): print "\n  ** Writing source data snippet to file: ** \n ", filename
	out_data.write(filename,format='fits')
	return filename

def to_catalogue_table(filename):
	"""
	Writes Aegean table to catalogue format (.csv) for plotting
	
	<param: filename> - the filename of the catalogue
	
	<return: N/A>
	"""
	
	# now redundant -> STILTS creates .csv table
	
	catalogue_filename = filename.replace(".fits","_catalogue.csv")
	if (verbose): print "\n <Writing catalogue to '.../{0[0]}/{0[1]}'>".format(catalogue_filename.split("\\")[-2:])
	
	data = Table.read(filename)
	data.write(catalogue_filename,format='ascii.csv')
	
def run_Aegean(fits_name, wide_catalogue_file, ra, dec, ang_diam, freq, path, Aegean_path):
	"""
	Runs Aegean.py source finding program by Paul Hancock. Outputs tables associated with the priorized source finding of the input .fits image
	
	<param: fits_name> - the input .fits filename to have sources detected with
	<param: wide_catalogue_file> - the input table of sources to be used in the priorized fitting
	<param: ra> - right ascension of the image
	<param: dec> - declination of the image
	<param: ang_diam> - angular diameter of the image
	<param: freq> - central frequency of the data
	<param: path> - the path to the input source table
	<param: Aegean_path> - the path to the directory with Aegean programs within
	
	<return: None>
	"""
	if (verbose): print "\n <Running Aegean.py>\n"
	
	out_filename = "GLEAM_catalogue_{0}_{1}_{2}_{3}".format(ra,dec,ang_diam,freq)
	input_filename = wide_catalogue_file.replace(".fits","")
	if (verbose): print "python \"{0}\\Aegean.py\" --input=\"{1}\" --priorized=1  --table=\"{2}\\{3}.fits\",\"{2}\\{3}.reg\" --telescope=MWA \"{4}\" ".format(Aegean_path,wide_catalogue_file,path,out_filename,fits_name)
	os.system("python \"{0}\\Aegean.py\" --input=\"{1}\" --priorized=1  --table=\"{2}\\{3}.fits\",\"{2}\\{3}.reg\" --telescope=MWA \"{4}\" ".format(Aegean_path,wide_catalogue_file,path,out_filename,fits_name))
	#os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\Aegean.py'+'"' + ' --input='+'"'+input_table_name+'"' + ' --priorized=1' + ' --table='+'"'+head+"\\"+out_filename+'.fits'+'","'+head+"\\"+out_filename+'.reg'+'"' + ' --telescope=MWA ' + '"'+input_fits_name+'"')
	
def run_AeRes(fits_file, catalogue_file, path, base, Aegean_path):
	"""
	Runs AeRes.py source subtracting program by Paul Hancock. Outputs a .fits image of the original image with the specified sources subtracted
	
	<param: fits_file> - the fits file to have sources subtracted from
	<param: catalogue_file> - the source catalogue of relevant extracted sources in Aegean format
	<param: path> - the path to the input source table
	<param: base> - the base name of the file #not implemented
	<param: Aegean_path> - the path to the directory with Aegean programs within
	
	<return: None>
	"""
	if (verbose): print "\n <Running AeRes.py>\n"
	
	
	
	if (verbose): print "python \"{0}/Aeres.py\" -c \"{1}\" -f \"{2}\" -r \"{3}/GLEAM_residual_{4}.fits\"".format(Aegean_path,catalogue_file,fits_file,path,base)
	os.system("python \"{0}/Aeres.py\" -c \"{1}\" -f \"{2}\" -r \"{3}/GLEAM_residual_{4}.fits\"".format(Aegean_path,catalogue_file,fits_file,path,base))
	#os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\AeRes.py'+'"' + ' -c ' + '"'+catalogue_filename+'"' + ' -f ' + '"'+fits_filename+'"' + ' -r ' + '"'+path+'\\GLEAM_residual_'+c_freq+'_'+base+'.fits"')

def run_BANE(fits_filename,Aegean_path):
	"""
	Runs Bane.py background and rms generator program by Paul Hancock. Outputs {bkg,rms}.fits files from the input .fits image.
	
	<param: fits_filename> - the file name of the fits file to have sources subtracted from
	<param: Aegean_path> - the path to the directory with Aegean programs within
	
	<return: None>
	"""
	if (verbose): print "\n <Running BANE.py>"
	
	print "python \"{0}\\BANE.py\" \"{1}\"".format(Aegean_path,fits_filename)
	os.system("python \"{0}\\BANE.py\" \"{1}\"".format(Aegean_path,fits_filename))
	#os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\BANE.py'+'"' + ' ' + '"'+fits_file+'"')
	
def get_cutout(ra, dec, freq, size=4.0, download_dir=None, listf=False):
	"""
	Automatically download GLEAM images from the postage stamp server using the template code that Chen has written.
	This function was written in majority by Paul Hancock, Aug-2015.
    
	<param: ra> - the centre RA of the map
	<param: dec> - the centre DEC of the map
	<param: freq> - central frequency of map; usage in file rename  
	<param: size> - the angular diameter of the map
	<param: download_dir> - Directory for which to save .fits image to
	<param: listf> - True/False depending on whether one wishes to print frequency list or not.
	
	<return: filename> - the file name of the downloaded .fitsfile
	"""
	if (verbose): print "\n <Downloading .fits file>"
	
	freq_ref = {'076':'072-080','084':'080-088','092':'088-095','099':'095-103','107':'103-111','115':'111-118','122':'118-126','130':'126-134','143':'139-147','151':'147-154','158':'154-162','166':'162-170','174':'170-177','181':'177-185','189':'185-193','197':'193-200','204':'200-208','212':'208-216','220':'216-223','227':'223-231','red':'072-103','green':'103-134','blue':'139-170','wide':'170-231'}
	try:
		freq_band = freq_ref[freq]
	except KeyError: # this should actually be handled by get_frequency()
		print " ** WARNING: no frequency '{0}' found **\n    Available frequencies: ".format(freq)
		for ii in freq_ref: print "      - {0}".format(ii)
		while True:
			choice = 'wide' if (auto) else str(raw_input("\n >> Choose frequency: ")) # auto_answer -> frequency = wide
			if (choice in freq_ref): freq_band = freq_ref[freq]; break
			else: print "\n  ** ERROR: invalid choice **  "
	
	gvp = GleamVoProxy() # start the gleam proxy // gvp = GleamVoProxy(p_port=7799)
	gvp.start()

	# check if downloads file exists, if not -> create it
	if (os.path.exists(download_dir) == False): 
		os.system("md \"{0}\"".format(download_dir))
		
	
	if (download_dir and (not os.path.exists(download_dir))):
		print "Invalid download dir: {0}".format(download_dir)
		return
	from pyvo.dal import sia
	svc = sia.SIAService(gvp.access_url) #start Simple Image Access service
	pos = (ra, dec) # position
	images = svc.search(pos, size)
	
	if listf:
		print "Available freq ranges are:"
		for img in images:
			print img.get('freq')
		return
	for img in images:
		# for each mached image, download or print its frequency and access url
		freq = img.get('freq')
		# only process the frequencies of interest
		if not freq in freq_band:
			continue
		print ' ** Downloading **'
		url = img.acref
		if (download_dir):
			download_file(url, ra, dec, freq, download_dir)
		else:
			print freq, url
	
	if (verbose): print "rename \"{0}\\{1}_{2}_{3}.fits\" \"GLEAM_cutout_{1}_{2}_{4}_{3}.fits\"".format(download_dir,ra,dec,freq_band,size)
	os.system("rename \"{0}\\{1}_{2}_{3}.fits\" \"GLEAM_cutout_{1}_{2}_{4}_{3}.fits\"".format(download_dir,ra,dec,freq_band,size))
	return "{0}\\GLEAM_cutout_{1}_{2}_{4}_{3}.fits".format(download_dir,ra,dec,freq_band,size)
	
	gvp.stop()
	
	
def main():
	usage = "usage: %prog [options] "
	parser = OptionParser(usage=usage)
	parser.add_option('-n', '--fits_file', 
					  action='store',type='string',dest='fits_file',default=None,
					  help=".fits file path", metavar="FITS_FILE")
	parser.add_option('-g','--galaxy',
					  action='store', type='string', dest='galaxy_name',default=None,
					  help="The name of the Dwarf galaxy",metavar="GALAXY_NAME")
	parser.add_option('-s','--source_data',
					  action='store', dest='data_file', default=None,
					  help="destination of input table for source data",metavar="SOURCE_DATA")
	parser.add_option('-f', '--freq', 
					  action='store',type='string',dest='freq', default="wide",
					  help="provide central frequency (MHz)", metavar="FREQUENCY")
	parser.add_option('-r','--ra',
					  action='store', type='float', dest='RA',default=None,
					  help="right ascension of the image",metavar="RA")
	parser.add_option('-d','--dec',
					  action='store', type='float', dest='DEC',default=None,
					  help="declination of the image",metavar="DEC")
	parser.add_option('-a','--angular_diameter',
					  action='store', type='float', dest='ang_diameter',default=4.0,
					  help="angular diameter of the sides of the image",metavar="ANGULAR_DIAMETER")
	parser.add_option('-w','--wighting',
					  action='store', type='string', dest='weighting',default=None,
					  help="Image weighting (-2 < R < +2)",metavar="WEIGHTING")
	parser.add_option('-b','--base',
					  action='store', type='string', dest='base_name',default=None,
					  help="The base name for the output Aegean formatted table")
	parser.add_option('-q', '--quiet',
					  action='store_false', dest='verbose', default=True,
					  help="don't print status messages to stdout")
	parser.add_option('-v','--verbose',
					  action='store_true', dest='verbose',default=False,
					  help="print status messages to stdout")
	parser.add_option('-y','--autoanswer',
					  action='store_true', dest='auto_answer',default=False,
					  help="auotmatically answer all choices")
	parser.add_option('-c','--clear_log',
					  action='store_true', dest='clear_log',default=False,
					  help="clear the log file")
	
	
	#parser.add_option("-","--",
	#				  action="", dest="",default=,
	#				  help="")
	
	print "\n  ********************************************************************\n  **************** REsidual Source IDentifier (RESID) ****************\n  ********************************************************************\n\n"
		
	(options, args) = parser.parse_args()
	global verbose; verbose = options.verbose
	if (verbose == True):
		print "  ** Verbose = True  **"
	else:
		print "  ** Verbose = False  **"
	global auto; auto = options.auto_answer
	if (verbose): print ("  ** Auto-Answer is on **  ")
	
	options.freq = get_frequency(options.freq)
	if (verbose): print "\n  ** Using frequency: {0}**".format(options.freq)

	# read log file, returns empty array if no log file exists or clear_log = True.
	if (options.clear_log == False): last_fits_file, last_gal_dir, last_gal_name, last_catalogue_file, last_Aegean_dir = read_log()
	else: os.system("rm ReSId_log.txt")
	
	if(verbose and options.clear_log != True):
		print "\n  ** Last used files: **"
		for kk in read_log(): print "    - {0}".format(kk) # print last usages
	
	#Tk().withdraw() # keeps then Tkinter root window from appearing
	
	# if no .fits file or galaxy has been specified
	if (options.galaxy_name == None and options.fits_file == None):
		if (last_fits_file != ''): print "\n Select an option:\n  1. Select a .fits file\n  2. Download .fits file from GLEAM server\n  3. Use previous .fits file: \"{0}\"".format(last_fits_file); num_choices = 3
		else: print "\n Select an option:\n  1. Select a .fits file\n  2. Download .fits file from GLEAM server\n"; num_choices = 2
		while True:
			try:
				choice = 3 if (auto) else int(raw_input(">>")) # auto_answer -> use previous file
				# choice = '3' if auto_answer==True else int(raw_input(">>"))
				if (choice>=1 and choice<=num_choices): break
				else: print "  ** ERROR: input out of bounds **"
			except ValueError:
				print "  ** ERROR: Invalid input **"
	
		if (choice == 1): # set fits_file to selected file
			while True:
				fits_file = askopenfilename(initialdir='\\'.join(last_fits_file.split('/')[0:-1])).replace('/','\\') # open dialog box and return the path to the selected file
				# AutoDownload option required
				if (fits_file != ''): break
				else: print "\n  ** ERROR: invalid selection **"
	
		if (choice == 2): # set fits_file to downloaded file
			if (options.RA == None or options.DEC == None): options.RA, options.DEC, options.ang_diameter = get_position(options.RA,options.DEC,options.ang_diameter)
			if (verbose): print " ** Downloading GLEAM cutout for \"{0}\" using parameters: **\n    - RA: {1}\n    - DEC: {2}\n    - Angular diameter: {3}\n    - Frequency: {4} MHz".format(galaxy,options.RA,options.DEC,options.ang_diameter,options.freq)
			if (os.path.exists("Downloads") == False): os.system("mkdir Downloads")
			################################# HERE ########################################
			
			os.mkdir("Downloads\\RA_{0}-DEC_{1}-FREQ_{2}-DIAM_{3}.fits".format(option.RA,options.DEC,options.freq,options.ang_diam))
			DL_file = get_cutout(options.RA, options.DEC, options.freq, options.ang_diameter, download_dir="Downloads\\RA_{0}-DEC_{1}-FREQ_{2}-DIAM_{3}.fits".format(option.RA,options.DEC,options.freq,options.ang_diam), listf=False)
			os.system("rename \"{0}\" \"GLEAM_cutout_{1}_{2}.fits\"".format(DL_file,freq,galaxy))
		
		if (choice == 3): # set fits_file to previous file
			fits_file = last_fits_file
	
	if (options.galaxy_name == None and options.fits_file == None): # *** Neither a 'fits_file', nor a 'galaxy_name' have been given
		if (last_fits_file != ''): # last_fits_file does exist
			print "\n ** WARNING: No .fits file specified **\n  previous .fits file used:\n \"{0}\"".format(last_fits_file)
			catch = False
			while True:
				if (catch): break
				choice = 'y' if (auto) else str(raw_input("\n>> Use this file (y/n)?:")) # auto_answer -> yes, use the previous fits_file
				if (choice.lower() == 'y'):
					fits_file = last_fits_file; catch = True # set fits_file to last selected 
				elif (choice.lower() == 'n'):
					print "  ** Select a GLEAM cutout .fits file ** "
					while True:
						fits_file = askopenfilename(initialdir='\\'.join(last_fits_file.split('/')[0:-1])).replace('/','\\') # open dialog box and return the path to the selected file
						# AutoDownload option required
						if (fits_file != ''): catch = True; break
						else: print "\n  ** ERROR: invalid selection **"
						
		else: # last_fits_file does not exist
			print " Select an option:\n  1. Download .fits file from GLEAM server\n  2. Select a .fits file"
			while True:
				try:
					choice = 1 if (auto) else int(raw_input(">>"))
					break
				except ValueError:
					print "  ** ERROR: Invalid input **"
			
			if (choice==1):
				if (options.RA == None or options.DEC == None): options.ra, options.dec, options.ang_diam = get_position(options.RA,options.DEC,options.ang_diameter)
				if (verbose): print " ** Downloading GLEAM cutout for \"{0}\" using parameters: **\n    - RA: {1}\n    - DEC: {2}\n    - Angular diameter: {3}\n    - Frequency: {4} MHz".format(galaxy,options.RA,options.DEC,options.ang_diameter,options.freq)
				DL_file = get_cutout(options.RA, options.DEC, options.freq, options.ang_diameter, download_dir=dir, listf=False)
				file = "{0}\\GLEAM_cutout_{1}_{2}.fits".format(dir,freq,galaxy)
				os.system("rename \"{0}\" \"GLEAM_cutout_{1}_{2}.fits\"".format(DL_file,freq,galaxy))
			elif (choice==2):
				print "  ** Select a GLEAM cutout .fits file ** "
				while True:	
					fits_file = askopenfilename().replace('/','\\') # open dialog box and return the path to the selected file
					if (fits_file != ''): break
					else: print "\n  ** ERROR: invalid selection **"
		
	# if both galaxy and .fits filename specified -> keep galaxy
	if (options.fits_file != None and options.galaxy_name != None): # *** 'fits_file' and 'galaxy_name' have been given.
		print "\n ** WARNING: '--fits_file' and '--galaxy' have been specified -> using '--galaxy' only **"
		options.fits_file = None
	
	# when galaxy name has been specified
	if (options.galaxy_name != None):
		last_gal_name = options.galaxy_name
		if (last_gal_dir != ""): # look in previous galaxy directory
			fits_file = find_gal_file(last_gal_dir,options.galaxy_name, options.RA, options.DEC, options.ang_diameter, options.freq)
			if (fits_file == ""):
				print "  ** WARNING: No GLEAM cutout was found for \"{0}\" **\n  ** Select a GLEAM cutout .fits file ** ".format(options.galaxy_name)
				while True:
					fits_file = askopenfilename(initialdir=last_gal_dir).replace('/','\\') # open dialog box and return the path to the selected file
					if (fits_file != ""): break
					else: print "\n  ** ERROR: invalid selection **"
				last_gal_dir = '\\'.join(fits_file.split('\\')[0:-2]) # update latest last_gal_dir
				last_gal_name = fits_file.split('\\')[-2]
		else:
			print "  ** Select the GLEAM cutout .fits file for the galaxy ** "
			fits_file = askopenfilename().replace('/','\\') # open dialog box and return the path to the selected file
			last_gal_dir = '\\'.join(fits_file.split('\\')[0:-2]) # update latest last_gal_dir
			last_gal_name = fits_file.split('\\')[-2]
	
	# if .fits filename has been specified
	if (options.fits_file != None):
		fits_file = find_file(options.fits_filename)
		if (fits_file == ""):
			print "  ** WARNING: No .fits file was found for \"{0}\" **\n  ** Select a GLEAM cutout .fits file ** ".format(options.fits_file)
			while True:
				fits_file = askopenfilename().replace('/','\\') # open dialog box and return the path to the selected file
				if (fits_file != ""): break
				else: print "\n  ** ERROR: invalid selection **"
	
	# specifying catalogue filename
	if (options.data_file != None): # if catalogue has been specified
		catalogue_file = options.data_file
	else: # if no catalogue has been specified
		if (last_catalogue_file != ""):
			print "\n ** WARNING: No catalogue file specified **\n  previous catalogue file used:\n \"{0}\" ".format(last_catalogue_file)
			catch = False
			while True:
				if (catch): break
				choice =  'y' if (auto) else str(raw_input("\n>> Use this file (y/n)?:"))
				if (choice.lower() == 'y'):
					catalogue_file = last_catalogue_file; catch = True # set catalgoue file to last used
				elif (choice.lower() == 'n'):
					print "  ** Select the input file for source data **  "
					while True:
						catalogue_file = askopenfilename(initialdir='\\'.join(last_catalogue_file.split('/')[0:-1])).replace('/','\\')
						if (catalogue_file != ""): catch = True; break
						else: print "  ** ERROR: invalid selection **"
		else:
			print "  ** Select the input file for source data **  "
			while True:
				catalogue_file = askopenfilename().replace('/','\\')
				if (catalogue_file != ""): break
				else: print "  ** ERROR: invalid selection **"
	
	
	# dir = path to directory, filename = name of the file only.
	dir, filename = '/'.join(fits_file.split('/')[:-1]), fits_file.split('/')[-1]
	print "ding ",dir.split("/")[-2:]
	if (verbose): print "\n  ** Using .fits file: '.../{0[0]}/{0[1]}/{1}' ** ".format(dir.split("/")[-2:],filename)
	
	if (verbose): print "\n  ** Getting .fits file header information  **"
	data, hdr = pyfits.getdata(fits_file,0,header=True) 
	Bmaj = hdr['BMAJ']*(60.0**2) # Beam major axis in arcsec
	Bmin = hdr['BMIN']*(60.0**2) # Beam minor axis in arcsec
	Bpa = hdr['BPA'] # Beam position angle
	# get robustness weighting
	if (options.weighting == None):
		Robust = hdr['ROBUST'] # Robustness of images. -1 = Briggs, 0 = Natural
	else
		Robust = options.weighting
	# get RA/DEC of center pixel if not already specified
	if (options.RA == None or options.DEC == None):
		if (verbose): print "  ** WARNING: No RA or DEC specified -> using RA and DEC from \".../{0}\"".format(filename)
		WCS = wcs.WCS(hdr, naxis=2)
		middle = WCS.wcs_pix2world(hdr['naxis1']/2,hdr['naxis2']/2,1) # get coordinates of middle pixel
		options.RA = middle[0]; options.DEC = middle[1]
		if (verbose): print "  RA = {0}\n  DEC = {1}".format(middle[0],middle[1])
	
	# make a base name for output files
	if (options.base_name == None):
		base = filename.replace(".fits","").replace("GLEAM_","").replace("Gleam_","").replace("cutout_","") if (options.galaxy_name==None) else "{0}_{1}".format(options.freq,options.galaxy_name)
	else:
		base = options.base_name
	
	# check for last Aegean directory
	if (last_Aegean_dir == ""):
		last_Aegean_dir = askdirectory(initialdir="~/",title="Select directory of 'Aegean' Tools").replace('/','\\')
		if (last_Aegean_dir == ""): print "\n  ** ERROR: AEGEAN directory was not given **\n   -- ABORTING --  "; exit() # exit if no directory given
		else: Aegean_dir = last_Aegean_dir
	else: Aegean_dir = last_Aegean_dir
	
	log(fits_file,last_gal_dir,last_gal_name,catalogue_file,Aegean_dir)
	
	# check if c_freq is RGB freq band -> in which case, need to generate table of found sources in .fits image by running Aegean
	if ((options.freq).lower() in ["red","r","green","g","blue","b"]):
		if (verbose): print "\n  ** Stacked frequency band: '{0}' chosen **  ".format(options.freq)
		run_BANE(fits_file,Aegean_dir) # *BANE.py does not work on windows ~Paul Hancock
		
		wide_catalogue_file = extract_sources(catalogue_file,options.RA,options.DEC,Bmaj,Bmin,Bpa,options.ang_diameter,"wide",dir,base)
		
		snippet_file = run_Aegean(fits_filename,wide_catalogue_file,options.RA,options.DEC,options.ang_diameter,options.freq,dir,Aegean_dir)
		os.system("rename \"{0}\\GLEAM_snippet_{1}_{2}_{3}_{4}_comp.fits\" \"GLEAM_snippet_{1}_{2}_{3}_{4}.fits\"".format(dir,options.RA,options.DEC,options.ang_diameter,options.freq))
		os.system("rename \"{0}\\GLEAM_snippet_{1}_{2}_{3}_{4}_comp.reg\" \"GLEAM_snippet_{1}_{2}_{3}_{4}.reg\"".format(dir,options.RA,options.DEC,options.ang_diameter,options.freq))
	else: # user has not specified an RGB frequency band	
		sources_file = extract_sources(catalogue_file,options.RA,options.DEC,Bmaj,Bmin,Bpa,options.ang_diameter,options.freq,dir,base)
			
	# run AeRes.py
	run_AeRes(fits_file, sources_file, dir, base, Aegean_dir)
	
	
	
main()
