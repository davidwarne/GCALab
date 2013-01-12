/* File: BitMapReader.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 18/04/2010
 *
 * Summary: Definitions of BitMap Reader Functions.
 *
 */


#ifndef BITMAPREADER_H
#define BITMAPREADER_H


#include "BitMapFile.h"

/*Core Reader Functions*/
ERROR BMP_ReadHeaders(BMPFILE*);
ERROR BMP_ReadColourPalette(BMPFILE*);
ERROR BMP_ReadImageData(BMPFILE*);

/*Higher order functions for reading, more user frendly*/
BMPImage ReadBMP(char*);

#endif
