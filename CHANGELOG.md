#changelog

# 2007-Jul-08 - Andrew W Stephens - alpha version
# 2007-Jul-09 - AWS - beta version
# 2007-Jul-10 - AWS - move most operations to cleanquad function
# 2007-Jul-11 - AWS - use stddev to decide whether to keep orig quad
# 2007-Jul-14 - AWS - generalized code to allow different pattern sizes
# 2007-Jul-18 - AWS - fix bug generating index arrays
# 2007-Jul-20 - AWS - add quadrant bias-level normalization
# 2007-Jul-23 - AWS - add force option
# 2007-Aug-06 - AWS - f/6 spectroscopy mode: use top & bottom for pattern
# 2007-Aug-22 - AWS - add -a flag to use all pixels
# 2008-Jan-11 - AWS - check for available image extensions
# 2008-Feb-05 - AWS - don't close input file until the end (req'd if next>2)
# 2008-Oct-02 - AWS - don't write pattern unless given the -p flag
# 2009-May-03 - AWS - use conformant default output file name
# 2009-May-13 - AWS - verify FITS header (old images have unquoted release date)
# 2009-May-22 - AWS - output full-frame pattern
# 2009-May-22 - AWS - improve quadrant bias normalization
# 2009-May-23 - AWS - add optional sky frame
# 2009-May-26 - AWS - add user-supplied bias offset correction
# 2009-Jul-04 - AWS - do not try to bias correct spectral flats
# 2009-Oct-24 - AWS - add basic row filtering
# 2009-Nov-06 - AWS - ignore bad pixels flagged in DQ extension
# 2009-Nov-08 - AWS - use mode for quadrant bias level normalization
# 2009-Nov-12 - AWS - use sigma-clipped stddev to judge quality of bias normalization
# 2009-Nov-17 - AWS - fit a Gaussian to the sky pixel distribution for bias norm.
# 2010-Feb-02 - AWS - sky subtract before quadrant normalization
# 2010-Feb-18 - AWS - add sigma-clipping to row filtering
# 2010-Apr-09 - AWS - only check for gcal status if OBSTYPE = FLAT
# 2010-Apr-10 - AWS - accept list input
# 2010-Apr-13 - AWS - minor tweak of the spectroscopic regions
# 2010-Jul-11 - AWS - allow images sizes which are multiples of the pattern size
# 2010-Oct-08 - AWS - option to read in bad pixel mask (e.g. object mask from nisky)
# 2010-Oct-10 - AWS - change the -g flag to take arguments
# 2010-Oct-11 - AWS - pad GNIRS images (2 row at the top)
# 2010-Oct-12 - AWS - GNIRS row filtering using an 8-pixel wide kernel
# 2010-Dec-21 - AWS - add grid filter
# 2010-Dec-28 - AWS - select GNIRS pattern region based on camera & slit
# 2011-Feb-03 - AWS - use extension 2 for nsprepared GNIRS data
# 2011-Feb-05 - AWS - add input glob expansion
# 2011-May-05 - AWS - output 32-bit files
# 2011-Jun-17 - AWS - catch modified GNIRS XD decker name
# 2012-Mar-29 - AWS - use the mode instead of median for pattern determination
# 2012-May-17 - AWS - do not modify input flag values
# 2013-Jun-16 - AWS - allow processing of non-standard FITS (with a warning)
# 2013-Jun-22 - AWS - add option to ignore DQ plane
# 2013-Jun-22 - AWS - use whole array if SKYIMAGE header keyword is present
# 2016-Jul-21 - AWS - propagate bad pixels through to the row filtering
# 2016-Oct-28 - AWS - switch pyfits -> astropy.io.fits
# 2017-Feb-01 - AWS - add option to reject pixels above a threshold before calculating pattern
# 2017-Mar-11 - AWS - fix calculation of stddev after quadrant adjustment
version = '2017-Mar-12' # astephens, code reorganization
version = '2018-Sep-09' # astephens, use numpy array masking, support fixing ROIs
version = '2018-Sep-11' # astephens, include quadrant leveling
version = '2018-Sep-13' # astephens, support manual user-supplied quadrant offsets
version = '2018-Sep-15' # astephens, playing with pattern shifts
version = '2018-Sep-23' # astephens, use multiprocessing
version = '2018-Sep-25' # astephens, -k to output the pattern kernel for testing
version = '2018-Oct-14' # astephens, use multiprocessing to level quadrants
version = '2018-Oct-19' # astephens, support user-supplied source ROIs to ignore
version = '2018-Oct-20' # astephens, add pshift parameter to change pattern shifts
version = '2018-Oct-21' # astephens, add row filtering
