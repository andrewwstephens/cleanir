#!/usr/bin/env python3

import argparse
from   astropy.io import fits
import datetime
import logging
import numpy
from   numpy import ma
import os
import sys

__version__ = '2021-Jun-13'

# --------------------------------------------------------------------------------------------------

def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.loglevel.upper()))
    logger = logging.getLogger()

    if len(args.fits) < 1:
        logger.error('No files found')
    else:
        for f in args.fits:
            row_subtract(f, args.roi, args.save)

# --------------------------------------------------------------------------------------------------

def row_subtract(fitsfile, rois, save):
    """
    Subtract the median of each row based on the pixels in the supplied ROIs.
    """
    logger = logging.getLogger('row_subtract')
    logger.info('Reading %s', fitsfile)
    hdulist = fits.open(fitsfile)
    ext = get_ext(hdulist)
    data = hdulist[ext].data
    mdata = ma.array(data, copy=True)
    ny, nx = data.shape
    logger.debug('Image size: %d x %d', nx, ny)
    roimask = ma.ones(data.shape)
    roimask.mask = True

    if rois is None:
        rois = ['1:%d' % nx]
    logger.debug('ROIs: %s', rois)

    for roi in rois:
        logger.debug('ROI: %s', roi)
        r = roi.split(':')
        if len(r) != 2:
            logger.error('ROIs must have 2 values: x1:x2:  %s', roi)
            raise SystemExit
        x1 = int(r[0]) - 1  # convert to zero-index
        x2 = int(r[1])      # zero index +1 because slicing does not include upper limit
        roimask.mask[0:ny, x1:x2] = False  # unmask the ROI: mask[y1:y2,x1:x2] = False
    logger.debug('roimask: %s', roimask)

    mdata *= roimask  # apply the ROI mask to the science data
    logger.debug('mdata: %s', mdata)
    median = ma.median(mdata, axis=1)  # calculate median of masked data along x-axis
    logger.debug('median: %s', median)
    logger.debug('len(median): %d', len(median))
    median_image = numpy.tile(median, nx).reshape(nx, ny).transpose()
    logger.debug('median_image: %s', median_image)
    logger.debug('median_image.shape: %s', median_image.shape)

    path, filename = os.path.split(fitsfile)
    if save is not None:
        for s in save:
            f = filename.replace('.fits', '_%s.fits' % s)
            logger.info('Writing %s', f)
            delete ([f])
            if s == 'masked':
                d = mdata.filled(fill_value=0)
            elif s == 'median':
                d = median_image.filled(fill_value=0)
            fits.PrimaryHDU(d).writeto(f)

    cleanfile = filename.replace('.fits', '_subtracted.fits')
    logger.info('Writing %s', cleanfile)
    delete([cleanfile])
    hdulist[ext].data = data - median_image.filled(fill_value=0)
    hdulist[0].header['ROWMED'] = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
    hdulist[0].header.comments['ROWMED'] = 'UT time stamp for rowmedian.py'
    hdulist.writeto(cleanfile, output_verify='silentfix')
    hdulist.close()
    return
    
# --------------------------------------------------------------------------------------------------

def get_ext(hdulist):
    logger = logging.getLogger('get_ext')
    numext = len(hdulist)
    logger.debug('Numext: %d', numext)
    sciext = None
    for i in range(len(hdulist)):
        try:
            extname = hdulist[i].header['EXTNAME']
        except:
            extname = None
        logger.debug('Extension %d name: %s', i, extname)
        if extname == 'SCI':
            sciext = i
            break
    if sciext is None:  # No extensions named 'SCI'
        if numext == 1:
            sciext = 0
        else:
            sciext = 1
    logger.debug('Science extension: %s', sciext)
    return sciext

# --------------------------------------------------------------------------------------------------

def delete(filelist):
    for f in filelist:
        if os.path.isfile(f):
            os.remove(f)
    return

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description='Subtract the median of each row calculated in the specified ROI(s).\n\n' +
        'examples:\n' +
        '  rowmedian N20210613S0123.fits --roi 300:400 --roi 600:700 --save median',
        epilog='Version: ' + __version__)
    parser.add_argument('fits', default=None, nargs='+', help='input FITS file(s)')
    parser.add_argument('--roi', default=None, type=str, action='append', metavar='X1:X2',
                        help='ROI(s) to use to calculate the row median.')
    parser.add_argument('--save', action='append', type=str, default=None,
                        choices=['masked', 'median'],
                        help='Save intermediate data product(s)')
    parser.add_argument('--loglevel', type=str, default='info',
                        choices=['debug', 'info', 'warning', 'error'], help='Log level [INFO]')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    main(args)

# --------------------------------------------------------------------------------------------------
