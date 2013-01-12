/* File: BitMapFile.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 19/04/2010
 *
 * Summary: Definitions of BitMap File Functions.
 *
 */
 
 
#ifndef BITMAPFILE_H
#define BITMAPFILE_H

#include<malloc.h>
#include "BitMapDefinitions.h"

BMPFILE* CreateBMPFILE(char*);
BMPFILE* CreateBMPFILE_FromImage(char*,BMPImage*);
ERROR DestroyBMPFILE(BMPFILE*);
ERROR BMP_OpenBitMap(BMPFILE*,char*);
ERROR BMP_CloseBitMap(BMPFILE*);

#endif
