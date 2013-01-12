/* File: BitMapWriter.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 19/04/2010
 *
 * Summary: Definitions of BitMap Writer Functions.
 *
 */
 
 
#ifndef BITMAPWRITER_H
#define BITMAPWRITER_H

#include "BitMapFile.h"

/*Core Writer functions*/
ERROR BMP_WriteHeaders(BMPFILE*);
ERROR BMP_WriteColourPalette(BMPFILE*);
ERROR BMP_WriteImageData(BMPFILE*);

/*higher level writer user friendly*/
ERROR WriteBMP(char*,BMPImage*);

#endif
