# CleanIR

CleanIR is an application to clean detector controller artifacts from [Gemini](https://www.gemini.edu) [NIRI](https://www.gemini.edu/sciops/instruments/niir) and [GNIRS](https://www.gemini.edu/sciops/instruments/gnirs) data.

## Table of Contents
- [Dependencies](#dependencies)
- [Features](#features)
- [Examples](#examples)
- [Support](#support)
- [Contributing](#contributing)

## Dependencies
CleanIR requires the following python packages:
- astropy - https://www.astropy.org/  
- matplotlib - http://matplotlib.org/  
- numpy - https://www.numpy.org/  
- scipy - https://www.scipy.org/  

## Features

CleanIR may be able to clean your IR data in several ways:

- **Full-frame Pattern Removal**
The full-frame "pattern" noise is the most common artifact affecting NIRI and GNIRS data.  It is noticeable by vertical 8-pixel wide stripes, or bright or dark columns with an 8-pixel periodicity.  This may be seen in one or more quadrants of the detector.  CleanIR analyzes each quadrant independently to calculate and subtract an 16x4 pixel noise kernel from each quadrant.

- **Partial-frame Pattern Removal**
A fraction of NIRI and GNIRS images are affected by the pattern noise over a section of the image.  In these cases the pattern will cover the full horizontal extent of one or more quadrants, but will either start or stop (or both) midway through the detector readout.  CleanIR accepts a command-line option to specify one or more regions that should be analyzed and cleaned.

- **Quadrant Leveling**
If your data have pattern noise there is a good chance that the bias level of the affected quadrants may not match the other quadrants.  CleanIR supports several methods of shifting the pattern noise as well as an option to match the background level of each quadrant.

- **Row Cleaning**
Less frequently NIRI and GNIRS data may be affected by row noise.  This is similar to the "pattern" noise in that it has an 8-pixel periodicity in a row, but it varies from row to row, unlike the normal "pattern" noise which affects all rows in a quadrant equally.  CleanIR can remove the component of this noise which is static across each row.

## Examples
  1. one
  2. two
  3. three

## Support
Feel free to contact astephens@gemini.edu with questions or problems.

## Contributing
Improvements are welcome!  Please follow the "fork-and-pull" Git workflow:
 1. **Fork** the repo on GitHub
 2. **Clone** the project to your own machine
 3. **Commit** changes to your own branch
 4. **Push** your work back up to your fork
 5. Submit a **Pull request**
