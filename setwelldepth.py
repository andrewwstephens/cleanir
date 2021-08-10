#!/usr/bin/env python3
import argparse
from astropy.io import fits
import sys
__version__ = '2021-Aug-10'  # astephens

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
    welldepth: the desired well depth, either None, 'shallow', or 'deep' (None prints the current values)
    """
    biasvoltage = {'deep': -0.87, 'shallow': -0.60}  # Detector bias voltage = VDDUC - VDET
    for i in expand(filerange):
        filename = 'N%sS%04d.fits' % (date,i)
        try:
            instrument = fits.getval(filename, 'INSTRUME')
            if instrument != 'NIRI':
                print(f'{filename} is {instrument}')
                continue
            vdduc1 = fits.getval(filename, 'VDDUC')
            vdduc2 = fits.getval(filename, 'A_VDDUC')
            vdet1  = fits.getval(filename, 'VDET')
            vdet2  = fits.getval(filename, 'A_VDET')
        except:
            print(f'{filename} Not found')
            continue
        bias1 = vdduc1 - vdet1
        bias2 = vdduc2 - vdet2

        if abs(bias1 - biasvoltage['shallow']) < 0.05 and abs(bias2 - biasvoltage['shallow']) < 0.05:
            current_welldepth = 'shallow'
        elif abs(bias1 - biasvoltage['deep']) < 0.05 and abs(bias2 - biasvoltage['deep']) < 0.05:
            current_welldepth = 'deep'
        else:
            current_welldepth = 'undefined'

        if welldepth is None:
            print(f'{filename} is {current_welldepth} well')
        elif welldepth == current_welldepth:
            print(f'{filename} is already {welldepth} well')
        else:
            newvdet1 = vdduc1 - biasvoltage[welldepth]
            newvdet2 = vdduc2 - biasvoltage[welldepth]
            print(f'{filename} is now {welldepth} well')
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
        description='Manually set the well depth in a NIRI FITS header.\n' +
                    'Note that you must be in the directory where the data are.\n\n' +
        'Examples:\n' + 
        'setwelldepth 20200822 1504-1509          -> print the current well depth\n' +
        'setwelldepth 20200822 1504-1509 deep     -> set these files to deep well\n' +
        'setwelldepth 20200822 1135-1138 shallow  -> set these files to shallow well\n',
        epilog='Version: ' + __version__)    
    parser.add_argument('date', type=str, help='Date of FITS files, e.g. 20201122')
    parser.add_argument('range', type=str, help='Range of file numbers, e.g. 1-3,5,7-9')
    parser.add_argument('welldepth', type=str, nargs='?', default=None,
                        choices=['shallow', 'deep'], help='Well depth')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    main(args)

# --------------------------------------------------------------------------------------------------
