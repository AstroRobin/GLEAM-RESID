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

# Imports for get_gleam()
from gleam_vo_example import GleamVoProxy, download_file
import pyvo

# Hardcoding
IDR_version = "4"

def choose(filenames):
	"""
	Given a list of filenames, this function will list all filenames in a formatted order and prompt the user to select a file
	
	<param: filenames> - a list of filenames
	
	<return: filename> - the chosen filename
	"""
	
	for kk in range(0,len(filenames)): print " {0} - {1}".format(kk+1,filenames[kk])
	while True:
			choice = raw_input("\n >> Choose file: ")
			try:
				choice = int(choice)
				if (choice >= 1 and choice <= len(filenames)):
					return filenames[choice-1]; break
				else:
					print "  ** ERROR: input out of bounds **  "
			except ValueError:
				print "  ** ERROR: invalid input **  "
	

def find_filename(search_path):
	"""
	Find .fits file in directory.
	
	:param search_path: The name of the fits file (and path) to be searched for
	
	:return: The .fits file name and path
	"""
	dir, filename = ntpath.split(search_path)
	if (verbose): print "\n <Searching for .fits file>\n  ** Searching for '{0}' in .../{1[0]}/{1[1]} **".format(filename,dir.split("\\")[-2:])
		
	found_filenames = []
	for file in os.listdir(dir):
		if filename in file:
			found_filenames.append(file)
	if (len(found_filenames) == 1):
		print "  ** Found .fits file: {0} **".format(found_filenames[0])
		filename = dir + "\\" + found_filenames[0]
	elif (len(found_filenames) > 1):
		print " ** Found multiple ({0}) .fits files ** ".format(len(found_filenames))
		filename = dir + "\\" + choose(found_filenames)
	else:
		print " ** No .fits files found with name '", filename,"' **\n    -- ABORTING --   "
		exit()
		
	return filename
	
def find_gal_filename(galaxy,ra,dec,ang_diam,freq):
	"""
	Find .fits file in directory for a dwarf galaxy
	
	:param filename: The name of the fits file to be searched for
	:param RA: right ascension of the map
	:param DEC: declination of the map
	:param ang_diam: angular diameter of the fits image
	:param freq: the frequency of the image
	
	:return: The .fits file name and path
	"""
	if (verbose): print "\n <Searching for .fits file>\n\n  ** Searching for galaxy: '{0}' **".format(galaxy)
	
	dir = "C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Dwarf Spheroidal Galaxies\\Images\\" # root directory for DSph galaxy images
	# look for directories with the name of the galaxy given
	catch = False
	gal_dirs = os.listdir(dir)
	for dir_name in gal_dirs: # iterate over all galaxy directories
		if (catch): break
		if (dir_name == galaxy):
			if (verbose): print "  ** Found directory: '{0}' **  \n".format(galaxy)
			dir = dir + dir_name
			catch = True; break
	if (catch == False):
		print "\n  ** WARNING: no directory '{0}' found **\n  ".format(galaxy)
		while True:
			choice = str(raw_input(" >> Make new galaxy directory '{0}' (y/n)?: ".format(galaxy)))
			if ("y" in choice.lower()): # make new directory for this folder
				os.system("mkdir \"{0}\"".format(dir+galaxy))
				print "  ** Directory: '{0}' has been created. **".format(galaxy)
				dir = dir + galaxy
				catch = True; break
			elif ("n" in choice.lower()): # don't make new directory -> ABORT
				if (verbose): print "\n    -- ABORTING --   "; exit()
	
	
	# look for files in galaxy directory with 'cutout' in their name
	if (verbose): print "  ** Searching for 'GLEAM_cutout' in '.../{0}' ** ".format(dir.split("\\")[-1])
	found_filenames = []
	gal_files = os.listdir(dir)
	for file_name in gal_files:
		if ("cutout" in file_name and freq in file_name):
			found_filenames.append(file_name)
	if (len(found_filenames) == 1): # if only one appropriate file found
		print "  ** Found .fits file: {0} **".format(found_filenames[0])
		filename = dir + "\\" + found_filenames[0]
	elif (len(found_filenames) > 1): # if multiple appropriate files found
		print "  ** Multiple ({0}) files found  ** ".format(len(found_filenames))
		filename = dir + "\\" + choose(found_filenames)
	else: # if no appropriate files found
		print " ** WARNING: No GLEAM_cutout_.fits files found in '{0}' directory **".format(galaxy)
		print " ** Download cutout for '{0}' using parameters: **\n    - RA: {1}\n    - DEC: {2}\n    - Angular diameter: {3}\n    - Frequency: {4}".format(galaxy,ra,dec,ang_diam,freq)
		while True:
			choice = "y" #str(raw_input(" >> Download (y/n)?: "))
			if ("y" in choice.lower()): # download cutout
				DL_filename = get_cutout(ra, dec, freq, ang_diam, download_dir=dir, listf=False)
				filename = "{0}\\GLEAM_cutout_{1}_{2}.fits".format(dir,freq,galaxy)
				os.system("rename \"{0}\" \"Gleam_cutout_{1}_{2}.fits\"".format(DL_filename,freq,galaxy))
				break
			elif ("n" in choice.lower()): # don't download cutout -> ABORT
				print "\n -- ABORTING --   "; exit()
			else: 
				print "\n  ** ERROR: invalid input **  "
				
	return filename

	
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
	
	if (verbose): print " <Finding frequency range>"	
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
	# function should check for existing filename outside of this function, i.e. in main()
	
	if (verbose): print '\n <Reading in data>'
	
	# this needs some cleaning up
	found_filename = check_for_file(head,RA,DEC,ang_diam)
	if (found_filename != None): filename = found_filename
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
	if (verbose): print "\n <Extracting sources> \n"
	
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
		if (data[ii]['RAJ2000'] >= RA - 0.5*ang_diam_buff and data[ii]['RAJ2000'] <= RA + 0.5*ang_diam_buff and data[ii]['DEJ2000'] >= DEC - 0.5*ang_diam_buff and data[ii]['DEJ2000'] <= DEC + 0.5*ang_diam_buff):
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
	# now redundant -> IDR4 has peak fluxes given.
	
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
	
	:return: the filename of the catalogue
	"""
	
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
	
	catalogue_filename = filename.replace(".fits","_catalogue.csv")
	if (verbose): print "\n <Writing catalogue to '.../{0[0]}/{0[1]}'>".format(catalogue_filename.split("\\")[-2:])
	
	data = Table.read(filename)
	data.write(catalogue_filename,format='ascii.csv')
	
	
def run_Aegean(input_fits_name, input_table_name, RA, DEC, ang_diam, freq, head):
	"""
	Runs Aegean.py source finding program by Paul Hancock. Outputs tables associated with the priorized source finding of the input .fits image
	
	:param input_fits_name: the input .fits filename to have sources detected with
	:param input_table_name: the input table of sources to be used in the priorized fitting
	:param RA: right ascension of the image
	:param DEC: declination of the image
	:param ang_diam: angular diameter of the image
	:param freq: central frequency of the data
	:param head: the path to the input source table
	"""
	if (verbose): print "\n <Running Aegean.py>"
	
	out_filename = 'GLEAM_snippet_'+RA+'_'+DEC+'_'+ang_diam+'_'+freq
	os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\Aegean.py'+'"' + ' --input='+'"'+input_table_name+'"' + ' --priorized=1' + ' --table='+'"'+head+"\\"+out_filename+'.fits'+'","'+head+"\\"+out_filename+'.reg'+'"' + ' --telescope=MWA ' + '"'+input_fits_name+'"')
	
def run_AeRes(fits_filename, catalogue_filename, c_freq, head, base):
	"""
	Runs AeRes.py source subtracting program by Paul Hancock. Outputs a .fits image of the original image with the specified sources subtracted
	
	:param fits_filename: the file name of the fits file to have sources subtracted from
	:param catalogue_filename: the file name of the source catalogue snippet in Aegean format
	:param c_freq: the central frequency.
	:param head: the path to the input source table
	:param base: the base name of the file #not implemented
	"""
	
	if (verbose): print "\n <Running AeRes.py>"
	print 'python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\AeRes.py'+'"' + ' -c ' + '"'+catalogue_filename+'"' + ' -f ' + '"'+fits_filename+'"' + ' -r ' + '"'+head+'\\GLEAM_residual_'+c_freq+'_'+base+'.fits"'
	os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\AeRes.py'+'"' + ' -c ' + '"'+catalogue_filename+'"' + ' -f ' + '"'+fits_filename+'"' + ' -r ' + '"'+head+'\\GLEAM_residual_'+c_freq+'_'+base+'.fits"')

def run_BANE(fits_filename):
	"""
	Runs Bane.py background and rms generator program by Paul Hancock. Outputs {bkg,rms}.fits files from the input .fits image.
	
	:param fits_filename: the file name of the fits file to have sources subtracted from
	"""
	
	if (verbose): print "\n <Running BANE.py>"
	os.system('python ' + '"'+'C:\\Users\\user\\OneDrive\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\Aegean\\Aegean-master\\BANE.py'+'"' + ' ' + '"'+fits_filename+'"')
	
def get_cutout(ra, dec, central_freq, size=4.0, download_dir=None, listf=False):
	"""
	Automatically download GLEAM images from the postage stamp server using the template code that Chen has written.
	This function was written in majority by Paul Hancock, Aug-2015.
    
	:param ra: the centre RA of the map
	:param dec: the centre DEC of the map
	:param central_freq: central frequency of map; usage in file rename  
	:param size: the angular diameter of the map
	:param download_dir: Directory for which to save .fits image to
	:param listf: True/False depending on whether one wishes to print frequency list or not.
	
	:return filename: the file name of the downloaded .fitsfile
	"""
	if (verbose): print "\n <Downloading .fits file>"
	
	freq_ref = {'076':'072-080','084':'080-088','092':'088-095','099':'095-103','107':'103-111','115':'111-118','122':'118-126','130':'126-134','143':'139-147','151':'147-154','158':'154-162','166':'162-170','174':'170-177','181':'177-185','189':'185-193','197':'193-200','204':'200-208','212':'208-216','220':'216-223','227':'223-231','red':'072-103','green':'103-134','blue':'139-170','wide':'170-231'}
	try:
		freqs = freq_ref[central_freq]
	except KeyError: # this should actually be handled by get_frequency()
		print " ** WARNING: no frequency '{0}' found **\n    Available frequencies: ".format(central_freq)
		for ii in freq_ref: print "      - {0}".format(ii)
		while True:
			choice = str(raw_input("\n >> Choose frequency: "))
			if (choice in freq_ref): freqs = freq_ref[central_freq]; break
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
		if not freq in freqs:
			continue
		print ' ** Downloading **'
		url = img.acref
		if (download_dir):
			download_file(url, ra, dec, freq, download_dir)
		else:
			print freq, url
	
	os.system('rename '+'"'+download_dir+"\\"+str(ra)+'_'+str(dec)+'_'+freqs+'.fits"'+' "GLEAM_cutout_'+str(ra)+'_'+str(dec)+'_'+str(size)+'_'+central_freq+'.fits"')
	return download_dir+'\\GLEAM_cutout_'+str(ra)+'_'+str(dec)+'_'+str(size)+'_'+central_freq+'.fits'
	
	gvp.stop()
	
	
def main():
	usage = "usage: %prog [options] "
	parser = OptionParser(usage=usage)
	parser.add_option('-n', '--fits_filename', 
					  action='store',type='string',dest='fits_filename',default=None,
					  help=".fits file name", metavar="FITS_FILE")
	parser.add_option('-g','--galaxy',
					  action='store', type='string', dest='galaxy_name',default=None,
					  help="The name of the Dwarf galaxy",metavar="GALAXY_NAME")
	parser.add_option('-q', '--quiet',
					  action='store_false', dest='verbose', default=True,
					  help="don't print status messages to stdout")
	parser.add_option('-v','--verbose',
					  action='store_true', dest='verbose',default=False,
					  help="print status messages to stdout")
	parser.add_option('-c','--catalogue',
					  action='store_true', dest='catalogue',default=False,
					  help="write to catalogue file")
	parser.add_option('-i','--datafile',
					  action='store', dest='data_filename', 
					  default="C:\\Users\\user\\OneDrive\\Documents\\Uni\\2016 - Semester 1\\Physics Dissertation\\GLEAM\\Data\\IDR{0}\\GLEAMIDR{0}.csv".format(IDR_version),
					  help="destination of input table for sources",metavar="SOURCES_FILE")
	parser.add_option('-f', '--central_freq', 
					  action='store',type='string',dest='central_freq', default="deep",
					  help="provide central frequency (MHz)", metavar="FREQUENCY")
	parser.add_option('-r','--ra',
					  action='store', type='float', dest='ra_map',default=None,
					  help="right ascension of the image",metavar="RA")
	parser.add_option('-d','--dec',
					  action='store', type='float', dest='dec_map',default=None,
					  help="declination of the image",metavar="DEC")
	parser.add_option('-a','--angular_diameter',
					  action='store', type='float', dest='ang_diameter',default=2.0,
					  help="angular diameter of the sides of the image",metavar="ANGULAR_DIAMETER")
	parser.add_option('-b','--base',
					  action='store', type='string', dest='base_name',default=None,
					  help="The base name for the output Aegean formatted table")
	
	
	#parser.add_option("-","--",
	#				  action="", dest="",default=,
	#				  help="")
	
	
	
	(options, args) = parser.parse_args()	
	global verbose; verbose = options.verbose
	
	options.central_freq = get_frequency(options.central_freq)
	if (verbose): print "\n  ** Using frequency: "+options.central_freq," **  "

	if (options.galaxy_name != None and options.fits_filename != None): 
		print "\n ** WARNING: Both -g (--galaxy) and -f (--fitsfile) have been specified **\n   -- ABORTING --   "
		exit()
	
	if (options.galaxy_name != None): # Galaxy FITS_Filename has been specified by user
		fits_filename = find_gal_filename(options.galaxy_name, options.ra_map, options.dec_map, options.ang_diameter, options.central_freq)
	else:
		if (options.fits_filename != None): # FITS_Filename has been specified by user
			fits_filename = find_filename(options.fits_filename)
		else:
			print "\n  ** No .fits filename specified **"
			if (options.ra_map == None or options.dec_map == None): print "  ** WARNING: Must specify both RA and DEC **\n   -- ABORTING --   "; exit()
			print "    Do you want to attempt downloading .fits file from < GLEAM Postage Stamp Service > using parameters: \n  - RA: {0}\n  - DEC: {1}\n  - Angular Diameter: {2}\n  - Frequency: {3} MHz".format(options.ra_map, options.dec_map, options.ang_diameter, options.central_freq)
			while True:
				choice = str(raw_input(">> (y/n)?: "))
				if (choice.lower() == 'y'):
					(out_dir_head, out_dir_tail) = ntpath.split(options.data_filename)
					fits_filename = get_cutout(options.ra_map, options.dec_map, options.central_freq, options.ang_diameter, download_dir="{0}\\Downloads".format(out_dir_head), listf=False)
					break
				elif (choice.lower() == 'n'):
					print "\n  ** Not attempting to download .fits filename **\n\n   -- ABORTING --   "
					exit()
					
					
			# fits_filename = str(raw_input(' >> Search for .fits file: '))

		
	head, tail = ntpath.split(fits_filename)
	if (verbose): print "\n  ** Using .fits file: '.../{0[0]}/{0[1]}/{0[2]}/{1}' ** ".format(head.split("\\")[-3:],tail)
	
	if (options.base_name == None):
		base = tail.replace(".fits","") if (options.galaxy_name==None) else options.galaxy_name
	else:
		base = options.base_name
	
	# check if c_freq is RGB freq band -> in which case, need to generate table of found sources in .fits image by running Aegean
	if ((options.central_freq).lower() in ["red","r","green","g","blue","b"]):
		# run_BANE(fits_filename) *BANE.py does not work on windows ~Paul Hancock
		if (verbose): print "\n  ** Stacked frequency band: '{0}' chosen **  ".format(options.central_freq)
		deep_filename = check_for_file(head,float(options.ra_map),float(options.dec_map),float(options.ang_diameter),"deep")
		
		if (deep_filename == None):
			# in future, do not exit -> rather automatically run through ReSId process and generate _deep snippet
			print "\n  ** Warning: No _deep GLEAM_snippet_...fits exists for this field **  \n  ** Please run ReSId.py with -central_freq=deep **\n   -- ABORTING --  "
			exit()
		
		run_Aegean(fits_filename,deep_filename,str(options.ra_map),str(options.dec_map),str(options.ang_diameter),options.central_freq,head)
		os.system("rename \"{0}\\GLEAM_snippet_{1}_{2}_{3}_{4}_comp.fits\" \"GLEAM_snippet_{1}_{2}_{3}_{4}.fits\"".format(head,options.ra_map,options.dec_map,options.ang_diameter,options.central_freq))
		os.system("rename \"{0}\\GLEAM_snippet_{1}_{2}_{3}_{4}_comp.reg\" \"GLEAM_snippet_{1}_{2}_{3}_{4}.reg\"".format(head,options.ra_map,options.dec_map,options.ang_diameter,options.central_freq))
	else: # user has not specified an RGB frequency band	
		# read data from input table
		in_data = read_data(options.data_filename,float(options.ra_map),float(options.dec_map),float(options.ang_diameter),head)
	
	catalogue_filename = check_for_file(head,float(options.ra_map),float(options.dec_map),float(options.ang_diameter),options.central_freq)
	if (catalogue_filename == None): # No '_snippet_' file was found
		# Extract sources which are constrained by input RA/DEC and ang_diam
		source_data = extract_sources(in_data,options.ra_map,options.dec_map,options.ang_diameter,head)
		
		# Convert source data to Aegean format table
		catalogue_filename = to_Aegean_table(source_data,options.central_freq,options.ra_map,options.dec_map,options.ang_diameter,head)
	if (options.catalogue): to_catalogue_table(catalogue_filename)
		
	#catalogue_filename = to_Aegean_table(source_data,options.central_freq,options.ra_map,options.dec_map,options.ang_diameter,head)
	
	# run AeRes.py
	run_AeRes(fits_filename, catalogue_filename, options.central_freq, head, base)
	
main()
