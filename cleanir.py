#!/usr/bin/env python3

import argparse
from   astropy.io import fits
from   astropy import stats
import datetime
from   functools import partial
import glob
import logging
from   matplotlib import pyplot
import multiprocessing
import numpy
from   numpy import ma
import os
from   scipy import optimize
from   scipy.stats import norm
import sys

version = '2019-Feb-19'

# --------------------------------------------------------------------------------------------------

# ToDo:
# Support --output so that ndisplay works with the 'clean' option!
# Update help
# Switch to python3
# Simpler to pass the four quadrants as a list, e.g. [q0,q1,q2,q3] ?
# Do we really need to save the variables qxsize and qysize ?
# Add regression tests
# Whenever there is intermittent noise, it seems to always be in coherent bands across the detector.
# This indicates that the pattern is actually across the *entire* detector, but with some offset or
# scaling between the quadrants.  If we could use all 4 quadrants (at the same time) for  pattern
# determination that may give better results, especially for GNIRS spectroscopy where only the
# unilluminated edges are available.
# -> either median combine the quadrants or the pattern, or filter the whole detector at once.

# --------------------------------------------------------------------------------------------------

def main(args):
    
    if args.debug:
        logger = ConfigureLogging('DEBUG')
    else:
        logger = ConfigureLogging()

    inputlist = glob.glob(args.fits)
    if len(inputlist) < 1:
        logger.error('No files found')
    logger.debug('inputlist = %s', inputlist)
    filelist = expand(inputlist)

    for f in filelist:
        cleanir = CleanIR(args)
        cleanir.clean(f)

    return

# --------------------------------------------------------------------------------------------------

class CleanIR():

    def __init__(self, args):
        logger = logging.getLogger('init')
        self.clip = args.clip
        self.debug = args.debug
        self.force = args.force
        self.graph = args.graph
        self.manual_q0_offset = args.TL
        self.manual_q1_offset = args.TR
        self.manual_q2_offset = args.BL
        self.manual_q3_offset = args.BR
        self.pshift = args.pshift
        self.pxsize = 16
        self.pysize =  4
        self.quadlevel = args.quadlevel
        self.roi = args.roi
        self.roimirror = True # mirror ROIs in the bottom quads to the top quads (and vice-versa?)
        self.rowfilter = args.rowfilter
        self.clip_sigma = args.sigma
        self.save = args.save
        self.src = args.src
        self.sub = args.sub
        self.use_dq = not args.nodq

    # ----------------------------------------------------------------------------------------------

    def clean(self, fitsfile):
        logger = logging.getLogger('clean')
        self.fitsfile = fitsfile
        self.read()
        self.mask_rois()
        self.mask_dq()
        self.subtract()
        self.mask_sources()
        self.calculate_pattern()
        self.subtract_pattern()
        self.row_filter()
        self.level_quadrants()
        self.write_output()
        return

    # ----------------------------------------------------------------------------------------------

    def read(self):
        logger = logging.getLogger('read')
        logger.info('Reading %s', self.fitsfile)

        self.hdulist = fits.open(self.fitsfile)

        numext = len(self.hdulist)
        logger.debug('Numext: %d', numext)

        if numext == 1:
            self.sciext = 0
        else:
            self.sciext = 1
        logger.debug('Science extension: %s', self.sciext)

        self.data = self.hdulist[self.sciext].data

        self.instrument = self.hdulist[0].header['INSTRUME']
        logger.debug('Instrument: %s', self.instrument)

        self.naxis1 = self.hdulist[self.sciext].header['naxis1']
        self.naxis2 = self.hdulist[self.sciext].header['naxis2']
        logger.debug('Image size: %s x %s', self.naxis1, self.naxis2)

        if self.instrument == 'NIRI':
            self.config = self.hdulist[0].header['FPMASK']

        elif self.instrument == 'GNIRS':
            self.config = self.hdulist[0].header['CAMERA'] + self.hdulist[0].header['DECKER']

            logger.debug('Padding GNIRS y-axis by 2 rows') # ydim must be a multiple of 4
            self.data = numpy.append(self.data, numpy.zeros((2,self.naxis1)), axis=0)
            self.naxis2 += 2
            logger.debug('New image size: %s x %s', self.naxis1, self.naxis2)

        else:
            logger.error('Unsupported instrument: %s', self.instrument)
            raise SystemExit

        logger.info('Config: %s', self.config)

        self.qxsize = int(self.naxis1 / 2)  # quadrant x size
        self.qysize = int(self.naxis2 / 2)  # quadrant y size

        self.mdata = ma.array(self.data, copy=True) # masked science data

        if self.instrument == 'GNIRS': # mask the padding
            self.mdata[-2:,] = ma.masked

        return

    # ----------------------------------------------------------------------------------------------

    def mask_rois(self):
        """
        Mask all pixels that are NOT in the user-supplied ROIs.
        These masked pixels will not be used in the pattern determination or be pattern subtracted.
        Expect ROI strings will have the format y1:y2
        If roimirror == True then apply ROIs to the opposite (top/bottom) half of the detector.
        """
        logger = logging.getLogger('mask_rois')
        if self.roi:
            logger.info('Masking around ROIs: %s', self.roi)
            self.roimask = ma.ones((self.naxis2, self.naxis1))
            self.roimask.mask = True
            for roi in self.roi:
                r = roi.split(':')
                if len(r) != 2:
                    logger.error('ROI must have 2 values: y1:y2:  %s', roi)
                    raise SystemExit
                y1 = int(r[0]) - 1 # convert to zero-index
                y2 = int(r[1])     # zero index +1 because slicing does not include upper limit
                logger.debug('...%d-%d', y1,y2)

                # Unmask the ROI: mask[y1:y2,x1:x2] = False
                self.roimask.mask[y1:y2,] = False

                if self.roimirror:
                    y3 = 1024-y2
                    y4 = 1024-y1
                    logger.debug('...%d-%d', y3,y4)
                    self.roimask.mask[y3:y4,] = False
 
            # Apply the ROI mask to the science data:
            self.mdata *= self.roimask

        return

    # ----------------------------------------------------------------------------------------------

    def mask_dq(self):
        """Mask pixels flagged in the DQ extension."""
        logger = logging.getLogger('mask_dq')
        
        if self.use_dq:
            try:
                dq = self.hdulist['DQ'].data
            except:
                dq = None
                logger.debug('No DQ extension found')

        else:
            dq = None

        if dq is not None:
            numpix = len(dq[dq>0])
            if numpix > 0:
                logger.info('Masking %d pixels flagged in DQ extension', numpix)
                self.dqmask = ma.masked_where(dq > 0, numpy.ones((self.naxis2, self.naxis1)))
                self.mdata *= self.dqmask
            else:
                self.dqmask = 1.0

        else:
            self.dqmask = 1.0
            
        return

    # ----------------------------------------------------------------------------------------------

    def subtract(self):
        logger = logging.getLogger('subtract')
        if self.sub:
            logger.info('Subtracting %s', self.sub)
            hdulist = fits.open(self.sub)
            numext = len(hdulist)
            logger.debug('...numext: %d', numext)
            if numext == 1:
                data_ext = 0
            else:
                data_ext = 1
            self.sub = hdulist[data_ext].data

            if self.instrument == 'GNIRS':
                logger.debug('Padding frame to subtract')
                self.sub = numpy.append(self.sub, numpy.zeros((2,self.naxis1)), axis=0)

            self.mdata -= self.sub
        else:
            self.sub = 0.0
        
        return
        
    # ----------------------------------------------------------------------------------------------

    def mask_sources(self):
        """
        Mask sources in the image by:
        1. Predefined NIRI and GNIRS spectroscopic ROIs
        2. User-supplied ROIs with the format x1:x2,y1:y2
        3. Iterative sigma clipping
        Masked regions will be ignored during pattern determination but still cleaned
        """
        logger = logging.getLogger('mask_sources')
        
        self.srcmask = ma.ones((self.naxis2, self.naxis1))
        
        if self.src is None: # The default, which is interpreted as "use the default mask"
                            
            if self.instrument == 'NIRI' and self.config in \
                   ('f6-2pixBl_G5214', 'f6-4pixBl_G5215', 'f6-6pixBl_G5216', 'f6-2pix_G5211'):
                logger.info('Using Y<=270 and Y>=728 for pattern determination')
                self.src = ['1:1024,271:727']

            elif self.instrument == 'GNIRS' and self.sub is not None and \
                     'Short' in self.config and not 'XD' in self.config:
                logger.info('Using X<=160 and X>=864 for pattern determination')
                self.src = ['161:863,1:1024']

            else:
                self.src = None

            if self.instrument == 'NIRI' and self.config in \
                     ('f6-4pix_G5212', 'f6-6pix_G5213', 'f32-6pix_G5229', 'f32-9pix_G5230'):
                logger.warning('Sky lines may be altered by pattern removal')

        logger.debug('Source mask: %s', self.src)

        if self.src is not None:
            if self.src[0] != 'None': # "None" is interpreted as do NOT use the default mask
                logger.info('Masking source ROI(s)...')
                for roi in self.src: # 'X1:X2,Y1:Y2'
                    r = roi.split(',')
                    x = r[0].split(':')
                    y = r[1].split(':')
                    x1 = int(x[0]) - 1 # convert to zero-index
                    x2 = int(x[1])
                    y1 = int(y[0]) - 1 # convert to zero-index
                    y2 = int(y[1])
                    logger.debug('%s -> %s %s %s %s', roi, x1, x2, y1, y2)
                    self.srcmask[y1:y2,x1:x2] = ma.masked

        if self.clip:
            logger.info('Sigma-clipping to remove sources...')
            q0,q1,q2,q3 = disassemble(self.data - self.sub)
            pool = multiprocessing.Pool(processes=4)
            mapfunc = partial(stats.sigma_clip,sigma_lower=3,sigma_upper=self.clip_sigma,maxiters=None)
            p0,p1,p2,p3 = pool.map(mapfunc, [q0,q1,q2,q3])
            nosrc = assemble(p0,p1,p2,p3) # the image with all sources masked
            sigmask = ma.ones((self.naxis2, self.naxis1))
            sigmask.mask = nosrc.mask
            self.srcmask *= sigmask

        self.mdata *= self.srcmask
        
        return

   # -----------------------------------------------------------------------------------------------

    def calculate_pattern(self):
        """
        Calculate the pattern for all four quadrants.
        """
        logger = logging.getLogger('calculate_pattern')
        logger.info('Calculating patterns...')

        # Create arrays of indices which correspond to the pattern tiled over the image
        indx = numpy.tile(numpy.arange(0, self.qxsize, self.pxsize), int(self.qysize/self.pysize))
        indy = numpy.arange(0, self.qysize, self.pysize).repeat(self.qxsize/self.pxsize)
        logger.debug('...indx: %s', indx)
        logger.debug('...indy: %s', indy)

        if self.graph is None:
            graph = False
        else:
            graph = True if 'pattern' in self.graph else False
    
        mapfunc = partial(cpq, pxsize=self.pxsize, pysize=self.pysize, qxsize=self.qxsize,
                          qysize=self.qysize, indx=indx, indy=indy, graph=graph, pshift=self.pshift)
        q0,q1,q2,q3 = disassemble(self.mdata)
        pool = multiprocessing.Pool(processes=4)
        #p0,p1,p2,p3 = pool.map(mapfunc, [q0,q1,q2,q3])
        (k0,p0),(k1,p1),(k2,p2),(k3,p3) = pool.map(mapfunc, [q0,q1,q2,q3])

        self.kernel  = assemble(k0,k1,k2,k3)
        self.pattern = assemble(p0,p1,p2,p3)

        if self.roi: # Zero regions not in the user-supplied ROIs
            self.pattern *= self.roimask.filled(fill_value=0)

        return

    # ----------------------------------------------------------------------------------------------

    def subtract_pattern(self):
        """
        Check and subtract the pattern from the original science image.
        """
        logger = logging.getLogger('subtract_pattern')
        logger.info('Checking pattern subtraction...')
        if self.force:
            logger.info('Forcing pattern subtraction on all quadrants')
        else:
            #d0,d1,d2,d3 = disassemble(self.data) # raw data
            d0,d1,d2,d3 = disassemble(self.mdata - self.sub) # masked subtracted data
            p0,p1,p2,p3 = disassemble(self.pattern)
            data = [(d0,p0), (d1,p1), (d2,p2), (d3,p3)] # pool.map only supports one input list
            pool = multiprocessing.Pool(processes=4)
            q0,q1,q2,q3 = pool.map(checksub, data)
            self.pattern = assemble(q0,q1,q2,q3)

        self.cleaned = self.data - self.pattern
        return

    # ----------------------------------------------------------------------------------------------

    def level_quadrants(self):
        """
        When full-frame (i.e. no user-supplied ROIs) compare each quadrant with the others.

        If there are user-supplied ROIs the region *outside* the ROIs is good and defines the
        reference bias level, so shift the bad regions to that level.
        """
        logger = logging.getLogger('level_quadrants')

        if self.quadlevel:

            logger.info('Leveling quadrants...')
            logger.debug('Fitting pixel distribution of each quadrant...')
            # Use the masked data to avoid bad pixels and objects.
            # Pattern subtract to minimize the dispersion.

            if self.graph is None:
                graph = False
            else:
                graph = True if 'offsets' in self.graph else False
            
            q0,q1,q2,q3 = disassemble(self.mdata - self.pattern)
            data = [ma.compressed(q0), ma.compressed(q1), ma.compressed(q2), ma.compressed(q3)]
            pool = multiprocessing.Pool(processes=4)
            mapfunc = partial(fitpixdist, graph=graph)
            g0,g1,g2,g3 = pool.map(mapfunc, data)            
            logger.debug('Peaks: %5.1f %5.1f %5.1f %5.1f', g0, g1, g2, g3)

            if self.roi:

                # Create a "not-ROI" mask which is the inverse of the ROI mask:
                notroi = ma.masked_array(self.roimask.data, ~self.roimask.mask)

                # The "good" background level is defined by the entire array excluding
                # pixels in the ROI, the DQ extension, and pixels with sources (stars, etc):
                # good = self.data * goodmask * self.dqmask * self.srcmask - self.sub
                # good = self.data * notroi - self.sub
                good = stats.sigma_clip(self.data * notroi - self.sub, 
                                        sigma_lower=3, sigma_upper=3.0, maxiters=1)
                middle = fitpixdist(ma.compressed(good), graph=self.graph)
                logger.debug('Good level: %.2f', middle)

            else: # full-frame
                # Should this be the *original* median of the entire image,
                # or the median after fixing the pattern noise?
                # Using the median after fixing pattern noise will give different
                # results for the different pattern shifting methods (min, mean, max).
                middle = numpy.median([g0,g1,g2,g3])
                logger.debug('Median of all quadrants: %.1f', middle)

            logger.debug('Calculating offsets...')
            o0 = middle - g0
            o1 = middle - g1
            o2 = middle - g2
            o3 = middle - g3
            logger.info('Quadrant offsets:  %.1f  %.1f  %.1f  %.1f', o0, o1, o2, o3)

        else:
            o0 = o1 = o2 = o3 = 0.0

        offsets = offset(numpy.zeros(self.naxis2*self.naxis1).reshape(self.naxis2,self.naxis1),
                         o0 + self.manual_q0_offset, o1 + self.manual_q1_offset,
                         o2 + self.manual_q2_offset, o3 + self.manual_q3_offset)

        if self.roi: # only apply offsets to the user-supplied ROI:
            offsets *= self.roimask.filled(fill_value=0)

        logger.debug('Applying offsets...')
        self.cleaned += offsets
        self.pattern -= offsets # the pattern is meant to be subtracted, so apply negative offsets
        self.kernel = offset(self.kernel,
                             -o0 - self.manual_q0_offset, -o1 - self.manual_q1_offset,
                             -o2 - self.manual_q2_offset, -o3 - self.manual_q3_offset)
        
        return

    # ----------------------------------------------------------------------------------------------

    def row_filter(self):
        """
        Generate an independent 8x1 pixel pattern for each row in each quadrant.
        """        
        logger = logging.getLogger('row_filter')

        if not self.rowfilter:
            return

        logger.info('Row filtering...')

        # quadrant AND/OR row filtering?
        # I'll have to find some examples to test, but IMHO we should try to remove as much of the
        # pattern noise on the full quadrant before trying to row filter, as the row filtering
        # will be affected by sky lines.  Alternatively we do quad filtering OR row filtering.
        
        q0,q1,q2,q3 = disassemble(self.mdata - self.sub - self.pattern)
        pool = multiprocessing.Pool(processes=4)
        p0,p1,p2,p3 = pool.map(rf, [q0,q1,q2,q3])
        self.rowpattern = assemble(p0,p1,p2,p3)

        self.cleaned -= self.rowpattern

        return
    
    # ----------------------------------------------------------------------------------------------

    def write_output(self):
        logger = logging.getLogger('write_output')

        path, filename = os.path.split(self.fitsfile)
        if len(path) > 0: path += '/'
        fileroot = path + 'c' + self.fitsfile[:self.fitsfile.rfind('.fits')]
        logger.debug('fileroot: %s', fileroot)

        if self.save is not None:
            for s in self.save:
                f = fileroot + '_' + s + '.fits'
                logger.info('Writing %s', f)
                delete ([f])
                if s == 'kernel':
                    data = self.kernel
                elif s == 'masked':
                    data = self.mdata.filled(fill_value=0)
                elif s == 'pattern':
                    data = self.pattern
                elif s == 'rowpattern':
                    data = self.rowpattern
                fits.PrimaryHDU(data).writeto(f)

        cleanfile = path + 'c' + self.fitsfile
        logger.info('Writing %s', cleanfile)
        delete ([cleanfile])

        if self.instrument == 'GNIRS':
            logger.debug('Removing padding')
            self.cleaned = numpy.delete(self.cleaned, [self.naxis2-1,self.naxis2-2], axis=0)

        self.hdulist[self.sciext].data = self.cleaned
        self.hdulist[0].header['CLEANIR'] = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
        self.hdulist[0].header.comments['CLEANIR'] = 'UT time stamp for cleanir.py'
        self.hdulist.writeto(cleanfile, output_verify='warn')
        self.hdulist.close()
        return
    
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

def disassemble(image):
    logger = logging.getLogger('disassemble')

    """
    Split an image into quadrants:
    +-------+
    | 0 | 1 |
    +---+---+
    | 2 | 3 |
    +---+---+
    """

    xsize = int(image.shape[1])
    ysize = int(image.shape[0])

    q0 = image[int(ysize/2):ysize,                   0:int(xsize/2)]
    q1 = image[int(ysize/2):ysize,        int(xsize/2):xsize]
    q2 = image[           0:int(ysize/2),            0:int(xsize/2)]
    q3 = image[           0:int(ysize/2), int(xsize/2):xsize]

    return q0, q1, q2, q3

# ----------------------------------------------------------------------------------------------

def assemble(q0, q1, q2, q3):
    """
    Assemble quadrants into an image.
    """
    logger = logging.getLogger('assemble')
    logger.debug('Assembling image...')

    xsize = 2 * q0.shape[1]
    ysize = 2 * q0.shape[0]

    if ma.is_masked(q0):
        image = ma.zeros(ysize*xsize).reshape(ysize,xsize)
    else:
        image = numpy.zeros(ysize*xsize).reshape(ysize,xsize)

    image[int(ysize/2):ysize,                   0:int(xsize/2)] = q0
    image[int(ysize/2):ysize,        int(xsize/2):xsize]        = q1
    image[           0:int(ysize/2),            0:int(xsize/2)] = q2
    image[           0:int(ysize/2), int(xsize/2):xsize]        = q3

    return image

# ----------------------------------------------------------------------------------------------

def offset(image, o0, o1, o2, o3):
    """
    Apply offsets to quadrants of an image.
    """
    logger = logging.getLogger('offset')
    logger.debug('Offsetting quadrant levels...')

    xsize = image.shape[1]
    ysize = image.shape[0]

    image[int(ysize/2):ysize,                   0:int(xsize/2)] += o0
    image[int(ysize/2):ysize,        int(xsize/2):xsize]        += o1
    image[           0:int(ysize/2),            0:int(xsize/2)] += o2
    image[           0:int(ysize/2), int(xsize/2):xsize]        += o3

    return image

# --------------------------------------------------------------------------------------------------

def checksub(data_and_pattern):
    """
    Calculate the 3-sigma clipped standard deviation of a quadrant before and after pattern
    subtraction and decide if there is improvement
    """
    logger = logging.getLogger('checksub')

    # The data and pattern arrays have been combined for pool.map; split them out here:
    data    = data_and_pattern[0]
    pattern = data_and_pattern[1]

    ostd = numpy.std(stats.sigma_clip(data,         sigma_lower=3.0, sigma_upper=3.0, maxiters=None))
    cstd = numpy.std(stats.sigma_clip(data-pattern, sigma_lower=3.0, sigma_upper=3.0, maxiters=None))

    if ostd - cstd > 0.01:
        logger.info('stddev %6.2f -> %6.2f  Improvement!', ostd, cstd)
    else:
        logger.info('stddev %6.2f -> %6.2f', ostd, cstd)
        pattern *= 0 # reset to zeros

    return pattern        

# --------------------------------------------------------------------------------------------------

def cpq(quad, pxsize=None, pysize=None, qxsize=None, qysize=None, indx=None, indy=None,
        graph=False, pshift='mean'):
    """
    Calculate the pattern for the supplied quadrant.
    """
    logger = logging.getLogger('cpq')

    # Create blank pattern array:
    p = numpy.zeros(pysize*pxsize).reshape(pysize,pxsize)

    median = ma.median(quad)
    stddev = ma.std(quad)
    logger.debug('...median: %.3f +/- %.3f', median, stddev)

    if graph:
        binwidth = 1
        binmin = median - 5. * stddev
        binmax = median + 5. * stddev
        bins = numpy.arange( binmin, binmax, binwidth )
        bincenters = bins[1:bins.size] - binwidth/2.
        iplot = 0

    for iy in range(0, pysize):
        for ix in range(0, pxsize):
            data = quad[indy+iy, indx+ix]
            p[iy,ix] = ma.median(data)

            if graph:
                iplot += 1
                plot = pyplot.subplot(pysize, pxsize, iplot)
                hist,bins = numpy.histogram(data, bins=bins)
                pyplot.plot(bincenters, hist, linestyle='-', marker='', markersize=1)
                pyplot.axvline(x=p[iy,ix], ls='--', color='red')
                if ix != 0:
                    plot.set_yticklabels([])

    if graph:
        logger.debug('...graphing results...')
        pyplot.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,wspace=0.,hspace=0.2)
        pyplot.show()

    # Decide how to shift the pattern:
    # For situations where the S/N is very low it is usually best to shift the pattern so that
    # it's mean is zero, and then subtracting the pattern has no net effect.  However, other
    # times this will consistently over-subtract the pattern, leaving dark quadrants.
    # In these situations it works better to shift the pattern so that it has a maximum of zero.

    # What is the best way to characterize the amplitude of the pattern?
    logger.debug('Pattern Amplitude: %.2f', numpy.amax(p) - numpy.amin(p))
    logger.debug('Pattern Sigma:     %.2f', numpy.std(p))


    if pshift == 'mean':
        logger.debug('Shifting the pattern to have a mean of zero')
        p -= numpy.mean(p)
    elif pshift == 'min':
        logger.debug('Shifting the pattern to have a minimum of zero')
        p -= numpy.amin(p)
    elif pshift == 'max':
        logger.debug('Shifting the pattern to have a maximum of zero')
        p -= numpy.amax(p)
    else:
        raise SystemExit('Invalid pshift')

    quadpattern = numpy.tile(p, (int(qysize/pysize), int(qxsize/pxsize))) # Tile pattern over quadrant

    return p, quadpattern

# --------------------------------------------------------------------------------------------------

def rf(quad): # Row filter a quadrant
    logger = logging.getLogger('rf')
    xsize = quad.shape[1]
    ysize = quad.shape[0]
    indx = numpy.arange(0, xsize, 8)
    pattern = ma.zeros((ysize,xsize))
    for iy in range(0, ysize):
        p = ma.zeros(8)
        for ix in range(0, 8):
            p[ix] = ma.median(quad[iy,indx+ix])
        pattern[iy] = numpy.tile(p, int(xsize/8))
    logger.debug('Row pattern mean: %s', ma.mean(pattern))
    pattern -= ma.mean(pattern) # set the mean to zero
    return pattern.filled(fill_value=0) # set masked values to zero

# --------------------------------------------------------------------------------------------------

def expand(listoflists):
    """Expand a list of file lists and return a list of files"""
    logger = logging.getLogger('expand')
    filelist = []
    for l in listoflists:
        if l.endswith('.fits'):
            filelist.append(l)
        else:
            logger.debug('Trying to expand %s', l)
            try:
                with open(l) as f:
                    for line in f:
                        filelist.append(line.strip())
            except:
                logger.error('Can\'t read %s', l)
                raise SystemExit

    filelist.sort()
    logger.debug('filelist = %s', filelist)
    return filelist

# --------------------------------------------------------------------------------------------------

def fitpixdist(array, graph=False):
    """
    Fit the pixel distribution of the passed array and return the center and width.
    Input:
        array = numpy masked array
    Output:
        center

    NOTE:  This works with v.1.14.5 but silently fails with numpy versions 1.15.0 and 1.15.1.
    """
    logger = logging.getLogger('fitpixdist')
    #logger.debug('array = %s', array)
 
    mu, sigma = norm.fit(array)
    logger.debug('mu = %.2f  sigma = %.2f', mu, sigma)

    mincts = numpy.amin(array)
    maxcts = numpy.amax(array)
    logger.debug('mincts = %.2f  maxcts = %.2f', mincts, maxcts)

    #bins = numpy.linspace(mincts, maxcts, 100)
    #logger.debug('bins = %s', bins)

    hist,bins = numpy.histogram(array, bins=int(maxcts-mincts))
    #logger.debug('Bins: %s', bins)

    bincenters = bins[1:bins.size] - (bins[1] - bins[0])/2.
    #logger.debug('Bin Centers: %s', bincenters)
    
    mode = bins[ hist.argmax() ]
    logger.debug('Mode = %.2f', mode)
    peak = hist.max()
    logger.debug('Peak = %.2f', peak)

    #fitsigma = 1.0 # how much to fit around the peak (+/- this sigma)
    #mincts = mode - fitsigma * inputstddev
    #maxcts = mode + fitsigma * inputstddev
    #t = bincenters[ (bincenters>mincts) & (bincenters<maxcts) ]
    #data =    hist[ (bincenters>mincts) & (bincenters<maxcts) ]

    #p0 = [mode, sigma, peak]
    #logger.debug('Initial parameter guesses: %.3f %.3f %.3f', p0[0], p0[1], p0[2])

    p0 = [mode, sigma, peak, 0.01, 0.00]
    logger.debug('Initial parameters:  %.3f %.3f %.3f %.3f %.3f', p0[0], p0[1], p0[2], p0[3], p0[4])

    xdata = bincenters
    ydata = hist
    #p, pcov = curve_fit(gaussian, xdata, ydata, p0=p0)
    p, pcov = optimize.curve_fit(gaussian_plus_slope, xdata, ydata, p0=p0)
    #print 'p =', p
    
    #logger.debug('Best fit parameters: %.3f %.3f %.3f', p[0], p[1], p[2])
    logger.debug('Best fit parameters: %.3f %.3f %.3f %.3f %.3f', p[0], p[1], p[2], p[3], p[4])

    xfit = xdata
    #yfit = gaussian(xfit, p[0], p[1], p[2])
    yfit = gaussian_plus_slope(xfit, p[0], p[1], p[2], p[3], p[4])

    if graph:
        pyplot.figure()
        pyplot.plot(xfit,  yfit,  linestyle='-', marker='', color='green')
        pyplot.plot(xdata, ydata, linestyle='',  marker='.', markersize=3, color='red')
        pyplot.show()

    return p[0]

# --------------------------------------------------------------------------------------------------

def gaussian(t,p0,p1,p2):    # p[0] = mu   p[1] = sigma    p[2] = peak
    return(p2 * numpy.exp( -(t - p0)**2 / (2 * p1**2) ))

# --------------------------------------------------------------------------------------------------

def gaussian_plus_slope(t,p0,p1,p2,p3,p4):  # p0=mu,  p1=sigma,  p2=peak,  p3=slope,  p4=offset
    return p2 * numpy.exp(-(t - p0)**2 / (2 * p1**2)) + p3 * t + p4

# --------------------------------------------------------------------------------------------------

def delete(filelist):
    for f in filelist:
        if os.path.isfile(f):
            os.remove(f)
    return

# --------------------------------------------------------------------------------------------------

def ConfigureLogging(level='INFO'):
    """Set up a console logger"""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) # set minimum threshold level for logger
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    consoleloghandler = logging.StreamHandler()
    if level.upper() == 'DEBUG':
        consoleloghandler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(name)-20s %(levelname)-8s %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')
    elif level.upper() == 'INFO':
        consoleloghandler.setLevel(logging.INFO)
    elif level.upper() == 'WARNING':
        consoleloghandler.setLevel(logging.WARNING)
    elif level.upper() == 'ERROR':
        consoleloghandler.setLevel(logging.ERROR)
    elif level.upper() == 'CRITICAL':
        consoleloghandler.setLevel(logging.CRITICAL)
    else:
        print ('ERROR: Unknown log error level')
        consoleloghandler.setLevel(logging.INFO)
    consoleloghandler.setFormatter(formatter)
    logger.addHandler(consoleloghandler)
    return logger

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='This script assumes that the NIRI/GNIRS pattern noise can be represented by a\
        16x4 pixel additive pattern which is repeated over the entire quadrant. The pattern is\
        determined for each quadrant by taking the mode of the pixel value distribution at each\
        position in the pattern.  Once the pattern has been determined for a\
        quadrant it is replicated to cover the entire quadrant and subtracted, and the mean of the\
        pattern is added back to preserve flux. The standard deviation of all the pixels in the\
        quadrant is compared to that before the pattern subtraction, and if no reduction was\
        achieved the subtraction is undone.  The pattern subtraction may be forced via the -f flag.\
        This process is repeated for all four quadrants and the cleaned frame is written to\
        c<infile> (or the file specified with the -o flag).  The pattern derived for each quadrant\
        may be saved with the -p flag.\
        \n\
        Pattern noise is often accompanied by an offset in the bias values between the four\
        quadrants.  One may want to use the -q flag to try to remove this offset.  This attempts to\
        match the iteratively determined median value of each quadrant. This method works best with\
        sky subtraction (i.e. with the -s flag), and does not work well if there are large extended\
        objects in the frame.  By default the median is determined from the entire frame, although\
        the -c flag will only use a central portion of the image.  Note that the derived quadrant\
        offsets will be applied to the output pattern file.\
        \n\
        Removing the pattern from spectroscopy is more difficult because of many vertical sky\
        lines.  By default NIRI f/6 spectroscopy with the 2-pixel or blue slits (which do not fill\
        the detector) use the empty regions at the bottom (1-272) and top (720-1024) of the array\
        for measuring the pattern.  This is not possible for other modes of spectroscopy where the\
        spectrum fills the detector. For these modes it is best to do sky subtraction before pattern\
        removal.  The quickest method is to pass a sky frame (or an offset frame) via the -s flag.\
        The manual method is to generate and subtract the sky, determine and save the pattern via\
        the -p flag, then subtract the pattern from the original image.  One may use the -a flag to\
        force using all of the pixels for the pattern determination.  If the SKYIMAGE FITS header\
        keyword is present it is assumed that the sky has already been subtracted and all pixels\
        will be used for the pattern determination.\
        \n\
        Note that you may use glob expansion in infile, however, the entire string must then be\
        quoted or any pattern matching characters (*,?) must be escaped with a backslash.',
        epilog='Version: ' + version)

    # Add comment about the pshift parameter.
    # Use "min" if the stripes are positive, and "max" if the stripes are negative,
    # and use "mean" if you don't want to change mean level of the image.

    parser.add_argument('fits', help='Fits file(s) or a list of FITS files')

    #parser.add_argument('-b', '--bpm', action='store', type=str, default=None,
    #                    help='Specify a bad pixel mask (overrides DQ plane)')

    parser.add_argument('--clip', action='store_true', default=False,
                        help='Iteratively sigma clip sources')

    #parser.add_argument('--dir', action='store', type=str, default=None,
    #                    help='Specify an input data directory')

    parser.add_argument('--debug', action='store_true', default=False, help='Debugging output')

    parser.add_argument('-f', '--force', action='store_true', default=False,
                        help='Force cleaning even if the stddev does not decrease')

    parser.add_argument('--graph', action='append', type=str, default=None,# metavar='OUTPUT',
                        choices=['pattern', 'offsets'],
                        help='Graph output')

    parser.add_argument('--nodq', action='store_true', default=False,
                        help='Ignore the DQ plane [False]')

    #parser.add_argument('-o', '--output', action='store', type=str, default=None,
    #                    help='Specify the cleaned output file [c<inputfile>]')

    parser.add_argument('--pshift', action='store', type=str, default='mean',
                        choices=['min', 'mean', 'max'],
                        help='Pattern shift.  Use "min" if the pattern appears as bright stripes on\
                        a dark background, "max" if the pattern appears as dark stripes on a bright\
                        background, or "mean" to not alter the original level of the image.')

    parser.add_argument('-q', '--quadlevel', action='store_true', default=False,
                        help='Level quadrants with additive offsets')

    # Use abbreviations to specify a quadrant (TBLR) to apply manual offsets:
    parser.add_argument('--TL', action='store', type=float, default=0.0, metavar='OFFSET',
                        help='Top Left quadrant manual offset (ADU)')
    parser.add_argument('--TR', action='store', type=float, default=0.0, metavar='OFFSET',
                        help='Top Right quadrant manual offset')
    parser.add_argument('--BL', action='store', type=float, default=0.0, metavar='OFFSET',
                        help='Bottom Left quadrant manual offset')
    parser.add_argument('--BR', action='store', type=float, default=0.0, metavar='OFFSET',
                        help='Bottom Right quadrant manual offset')

    parser.add_argument('--roi', default=None, type=str, action='append', metavar='Y1:Y2',
                        help='ROI(s) to clean')

    parser.add_argument('--rowfilter', action='store_true', default=False,
                        help='Filter each quadrant row with an 8-pixel kernel.  This may be useful\
                        for GNIRS XD spetra that cannot be cleaned otherwise.')

    parser.add_argument('--save', action='append', type=str, default=None,
                        choices=['kernel', 'masked', 'pattern', 'rowpattern'],
                        help='Save intermediate product(s)')

    parser.add_argument('--sigma', action='store', type=float, default=3.0,
                        help='Source clipping upper limit [3.0]')

    parser.add_argument('--src', default=None, type=str, action='append', metavar='X1:X2,Y1:Y2',
                        help='ROI(s) to ignore when calculating the pattern.  Note that default\
                        ROIs will be used for some spectroscopic configurations.  Set to "None"\
                        to NOT use the default ROIs.')
    # --src defaults to None, which this script interprets as "use the default mask".
    # If the user wants to NOT use the default mask they should pass in the string "None"
    # which will come through as args.src = ['None'].
    
    parser.add_argument('--sub', action='store', type=str, default=None, metavar='FITS',
                        help='FITS file to subtract before calculating pattern')

    args = parser.parse_args()
    main(args)

# --------------------------------------------------------------------------------------------------
