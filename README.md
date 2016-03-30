# GLEAM-RESID
GLEAM REsidual Source IDentifier

This program is designed to IDentify REsidual Sources that were possibly missed in the initial GLEAM source identification selected from the WHITE frequency stacked (170-231MHz) images in IDR3. Such sources are believed to be very steep spectrum sources (alpha > 2) such that they could be detected at low frequencies images, yet fall below the detection threshold in the WHITE image frequency range. 

Options:
  -h, --help            show this help message and exit
  -f FILE, --file=FILE  write report to FILE
  -q, --quiet           don't print status messages to stdout
  -c FREQUENCY, --centralfreq=FREQUENCY
                        provide central frequency
  -i SOURCES_FILE, --datafile=SOURCES_FILE
                        destination of input table for sources
  -r RA, --ra=RA        right ascension of the image
  -d DEC, --dec=DEC     declination of the image
  -a ANGULAR_DIAMETER, --angulardiameter=ANGULAR_DIAMETER
                        angular diameter of the sides of the image
  -b BASE_NAME, --base=BASE_NAME
                        The base name for the output Aegean formatted table
