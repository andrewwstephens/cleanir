# changelog

2007-Jul-08, astephens, alpha version  
2007-Jul-10, astephens, move most operations to cleanquad function  
2007-Jul-11, astephens, use stddev to decide whether to keep orig quad  
2007-Jul-14, astephens, generalized code to allow different pattern sizes  
2007-Jul-18, astephens, fix bug generating index arrays  
2007-Jul-20, astephens, add quadrant bias-level normalization  
2007-Jul-23, astephens, add force option  
2007-Aug-06, astephens, f/6 spectroscopy mode: use top & bottom for pattern  
2007-Aug-22, astephens, add -a flag to use all pixels  
2008-Jan-11, astephens, check for available image extensions  
2008-Feb-05, astephens, don't close input file until the end (req'd if next>2)  
2008-Oct-02, astephens, don't write pattern unless given the -p flag  
2009-May-03, astephens, use conformant default output file name  
2009-May-13, astephens, verify FITS header (old images have unquoted release date)  
2009-May-22, astephens, output full-frame pattern  
2009-May-22, astephens, improve quadrant bias normalization  
2009-May-23, astephens, add optional sky frame  
2009-May-26, astephens, add user-supplied bias offset correction  
2009-Jul-04, astephens, do not try to bias correct spectral flats  
2009-Oct-24, astephens, add basic row filtering  
2009-Nov-06, astephens, ignore bad pixels flagged in DQ extension  
2009-Nov-08, astephens, use mode for quadrant bias level normalization  
2009-Nov-12, astephens, use sigma-clipped stddev to judge quality of bias normalization  
2009-Nov-17, astephens, fit a Gaussian to the sky pixel distribution for bias norm.  
2010-Feb-02, astephens, sky subtract before quadrant normalization  
2010-Feb-18, astephens, add sigma-clipping to row filtering  
2010-Apr-09, astephens, only check for gcal status if OBSTYPE = FLAT  
2010-Apr-10, astephens, accept list input  
2010-Apr-13, astephens, minor tweak of the spectroscopic regions  
2010-Jul-11, astephens, allow images sizes which are multiples of the pattern size  
2010-Oct-08, astephens, option to read in bad pixel mask (e.g. object mask from nisky)  
2010-Oct-10, astephens, change the -g flag to take arguments  
2010-Oct-11, astephens, pad GNIRS images (2 row at the top)  
2010-Oct-12, astephens, GNIRS row filtering using an 8-pixel wide kernel  
2010-Dec-21, astephens, add grid filter  
2010-Dec-28, astephens, select GNIRS pattern region based on camera & slit  
2011-Feb-03, astephens, use extension 2 for nsprepared GNIRS data  
2011-Feb-05, astephens, add input glob expansion  
2011-May-05, astephens, output 32-bit files  
2011-Jun-17, astephens, catch modified GNIRS XD decker name  
2012-Mar-29, astephens, use the mode instead of median for pattern determination  
2012-May-17, astephens, do not modify input flag values  
2013-Jun-16, astephens, allow processing of non-standard FITS (with a warning)  
2013-Jun-22, astephens, add option to ignore DQ plane  
2013-Jun-22, astephens, use whole array if SKYIMAGE header keyword is present  
2016-Jul-21, astephens, propagate bad pixels through to the row filtering  
2016-Oct-28, astephens, switch pyfits -> astropy.io.fits  
2017-Feb-01, astephens, add option to reject pixels above a threshold before calculating pattern  
2017-Mar-11, astephens, fix calculation of stddev after quadrant adjustment  

---
**Current code-base starts here.**  
2017-Mar-12, astephens, complete code reorganization  
2018-Sep-09, astephens, use numpy array masking to support fixing ROIs  
2018-Sep-11, astephens, include quadrant leveling  
2018-Sep-13, astephens, support manual user-supplied quadrant offsets  
2018-Sep-23, astephens, use multiprocessing  
2018-Sep-25, astephens, add -k option to output the pattern kernel for testing  
2018-Oct-14, astephens, use multiprocessing to level quadrants  
2018-Oct-19, astephens, support user-supplied source ROIs to ignore  
2018-Oct-20, astephens, add pshift parameter to change pattern shifts  
2018-Oct-21, astephens, add row filtering  
2019-Feb-10, astephens, update to python3
2019-Feb-19, astephens, fix missing int in row filter
2019-Apr-08, astephens, close and join pools when finished
