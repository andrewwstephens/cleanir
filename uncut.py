#!/usr/bin/env python3
import argparse
from   astropy.io import fits
import logging
import numpy
import os
import re
import sys
__version__ = '2021-May-29'  # astephens, beta version

# --------------------------------------------------------------------------------------------------

def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.loglevel.upper()))
    logger = logging.getLogger()
    if len(args.fits) < 1:
        logger.error('No files found')
    else:
        for f in args.fits:
            uncut(f)

# --------------------------------------------------------------------------------------------------

def uncut(fitsfile):
    """Undo NSCUT"""
    logger = logging.getLogger('uncut')
    logger.info('Reading %s', fitsfile)
    hdulist = fits.open(fitsfile)
    if 'NSCUT' not in hdulist[0].header:
        logger.info('image has not been NSCUT')
        return
    logger.info('Cut image size: %d x %d', hdulist['SCI'].data.shape[1], hdulist['SCI'].data.shape[0])
    uncutfilename = fitsfile.replace('.fits', '_uncut.fits')
    delete([uncutfilename])
    for i in range(1, len(hdulist)):  # the PHU at [0] does not have an XTENSION keyword
        xtension = hdulist[i].header['XTENSION']
        logger.debug('%d: %s', i, xtension, )
        if xtension == 'IMAGE':
            hdulist[i].data = pad(hdulist[i]).data
    logger.info('Uncut image size: %d x %d', hdulist['SCI'].data.shape[1], hdulist['SCI'].data.shape[0])
    logger.info('Writing %s', uncutfilename)
    hdulist.writeto(uncutfilename, output_verify='silentfix')  # {ignore, fix, silentfix, warn}
    hdulist.close()
    return

# --------------------------------------------------------------------------------------------------

def pad(image):
    """Pad an NSCUT GNIRS image back to it's original size."""
    logger = logging.getLogger('pad')
    nscutsec = image.header['NSCUTSEC']
    logger.debug('NSCUTSEC: %s', nscutsec)
    x1, x2, y1, y2 = map(int, re.split('[:,]', nscutsec.strip('[]')))
    # uncut GNIRS images are [1:1024,1:1022]
    pad_x_front = x1 - 1
    pad_x_end = 1024 - x2
    pad_y_front = y1 - 1
    pad_y_end = 1022 - y2
    logger.debug('X padding: %d, %d', pad_x_front, pad_x_end)
    logger.debug('Y padding: %d, %d', pad_y_front, pad_y_end)
    extname = image.header['EXTNAME']
    logger.debug('EXTNAME: %s', extname)
    if extname == 'DQ':
        padval = 1
    else:
        padval = 0
    logger.debug('Padding value: %d', padval)
    data = image.data
    if pad_y_front > 0:
        data = numpy.append(numpy.ones((pad_y_front, data.shape[1]), dtype=numpy.int8) * padval, data, axis=0)
    if pad_y_end > 0:
        data = numpy.append(data, numpy.ones((pad_y_end, data.shape[1]), dtype=numpy.int8) * padval, axis=0)
    if pad_x_front > 0:
        data = numpy.append(numpy.ones((data.shape[0], pad_x_front), dtype=numpy.int8) * padval, data, axis=1)
    if pad_x_end > 0:
        data = numpy.append(data, numpy.ones((data.shape[0], pad_x_end), dtype=numpy.int8) * padval, axis=1)
    return data

# --------------------------------------------------------------------------------------------------

def delete(filelist):
    for f in filelist:
        if os.path.isfile(f):
            os.remove(f)
    return

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Undo an NSCUT so that a GNIRS image may be processed by cleanir.py',
        epilog='Version: ' + __version__)
    parser.add_argument('fits', default=None, nargs='+', help='input FITS file(s)')
    parser.add_argument('--loglevel', type=str, default='info',
                        choices=['debug', 'info', 'warning', 'error'], help='Log level [INFO]')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    main(args)

# --------------------------------------------------------------------------------------------------
