/* File: BitMapDefintions.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 18/04/2010
 *
 * Summary: Definitions of BitMap Datafile Format, reader functions use the BMPFILE
 *          structure for all read and write operations (see BitMapReader.h and BitMapWriter.h)
 */

#ifndef BITMAPDEFINITIONS_H
#define BITMAPDEFINITIONS_H

/* Compression level macros*/
#define BI_RGB 0
#define BI_RLE8 1
#define BI_RLE4 2
#define BI_BITFIELDS 3
#define BI_JPEG 4
#define BI_PNG 5


typedef int ERROR;

#define NO_ERRORS 0
#define FILE_OPEN_ERROR 1
#define FILE_READ_ERROR 2
#define MEMORY_ERROR 3
#define FILE_CLOSE_ERROR 4
#define FILE_WRITE_ERROR 5

#include <stdint.h>
#include <stdio.h>

/*BitMap File Header*/
typedef struct 
{
	unsigned char type[2];
	uint32_t fileSize;
	uint16_t reserved1;
	uint16_t reserved2;
	uint32_t offset;
}BMPFILEHEADER;


/*The most common Windows 3.0 BitMap format info header*/
typedef struct
{
	uint32_t size;
	uint32_t width;
	uint32_t height;
	uint16_t planes;
	uint16_t colorDepth;
	uint32_t compression;
	uint32_t imageSize;
	uint32_t hRes;
	uint32_t wRes;
	uint32_t paletteSize;
	uint32_t numImportantColours;
	
}BMPINFOHEADER;

/* these are to handel support for newer BitMao formats*/
typedef struct
{
	
}BMPINFOHEADERV4;

typedef struct
{
	
}BMPINFOHEADERV5;

/* for use with Bitmaps with colour depth less than 24 bit*/
typedef struct
{
	int size;
	uint32_t* colours;
}PALETTE;

/* Bitmap file structure*/
typedef struct
{
	char* fileName;
	FILE* bmpfile_fp;
	BMPFILEHEADER* bmpFileHeader;
	BMPINFOHEADER* bmpInfoHeader;
	PALETTE* bmpPalette;
	uint8_t* imageData;
}BMPFILE;

typedef struct
{
	unsigned long int width;
	unsigned long int height;
	unsigned char** RGB; /*each row i column j -> RGB[i][j],RGB[i][j+1],RGB[i][j+2] */
}BMPImage;

#endif
