#!/usr/bin/env python3

import argparse
from   astropy.io import fits

__version__ = '2020-Nov-22'  # Andrew Stephens

# --------------------------------------------------------------------------------------------------

def main(args):
    setwelldepth(args.date, args.range, args.welldepth)
    return

# --------------------------------------------------------------------------------------------------

def setwelldepth(date, filerange, welldepth):
    """
    Manually set the well depth in a NIRI FITS header.
    date: date of the FITS file, e.g. 20201122
    filerange: range of file numbers, e.g. 1-3,5,7-9
    welldepth: the desired well depth, either 'shallow' or 'deep'
    """
    biasvoltage = {'deep': -0.87, 'shallow': -0.60}  # Detector bias voltage = VDDUC - VDET
    for i in expand(filerange):
        filename = 'N%sS%04d.fits' % (date,i)
        try:
            vdduc1 = fits.getval(filename, 'VDDUC')
            vdduc2 = fits.getval(filename, 'A_VDDUC')
            vdet1  = fits.getval(filename, 'VDET')
            vdet2  = fits.getval(filename, 'A_VDET')
        except:
            print('%s  Not found' % filename)
            continue
        bias1 = vdduc1 - vdet1
        bias2 = vdduc2 - vdet2
        if  abs(bias1 - biasvoltage[welldepth]) < 0.05 and \
            abs(bias2 - biasvoltage[welldepth]) < 0.05:
            #print('%s  %.3f - %.3f = %.3f %s' % (filename, vdduc1, vdet1, bias1, welldepth.upper()))
            print('%s is already set to %s well' % (filename, welldepth.upper()))
        else:
            newvdet1 = vdduc1 - biasvoltage[welldepth]
            newvdet2 = vdduc2 - biasvoltage[welldepth]
            #print('%s  %.3f - %.3f = %.3f --> setting VDET = %.2f' % (filename, vdduc1, vdet1, bias1, newvdet1))
            print('%s is now set to %s well' % (filename, welldepth.upper()))
            fits.setval(filename, 'VDET',   value=newvdet1)
            fits.setval(filename, 'A_VDET', value=newvdet2)
    return

# --------------------------------------------------------------------------------------------------

def expand(list_of_ranges):
    """Expand a list of ranges and return a python list, e.g. 1-3,5,7-9 -> [1,2,3,5,7,8,9]"""
    output = []
    for r in list_of_ranges.split(','):
        if '-' in r:
            indx = r.find('-')
            output.extend(list(range(int(r[:indx]),int(r[1+indx:])+1)))
        else:
            output.extend([int(r)])
    return output

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description='Manually set the well depth in a NIRI FITS header.\n\n' + 
        'Examples:\n' + 
        'setwelldepth 20200822 1504-1509 deep\n' +
        'setwelldepth 20200822 1135-1138 shallow\n',
        epilog='Version: ' + __version__)    
    parser.add_argument('date', type=str, help='Date of FITS files, e.g. 20201122')
    parser.add_argument('range', type=str, help='Range of file numbers, e.g. 1-3,5,7-9')
    parser.add_argument('welldepth', type=str,  choices=['shallow', 'deep'], help='well depth')
    args = parser.parse_args()
    main(args)

# --------------------------------------------------------------------------------------------------
