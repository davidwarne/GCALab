/* File: BitMapReader.cpp
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 18/04/2010
 *
 * Summary: Implementations of BitMap Reader Functions.
 *
 * NOTE: at this point on uncompressed BITMAPS Supported
 */
 
 
#include "BitMapReader.h"

ERROR BMP_ReadHeaders(BMPFILE* bmp_fp)
{
	bmp_fp->bmpFileHeader = (BMPFILEHEADER*)malloc(sizeof(BMPFILEHEADER));
	bmp_fp->bmpInfoHeader = (BMPINFOHEADER*)malloc(sizeof(BMPINFOHEADER));
	if(!(fread((void*)(bmp_fp->bmpFileHeader->type),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpFileHeader->fileSize),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpFileHeader->reserved1),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpFileHeader->reserved2),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpFileHeader->offset),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	
	fseek(bmp_fp->bmpfile_fp,0x000E,SEEK_SET);
	
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->size),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->width),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->height),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->planes),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->colorDepth),2,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->compression),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->imageSize),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->hRes),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->wRes),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->paletteSize),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	if(!(fread((void*)&(bmp_fp->bmpInfoHeader->numImportantColours),4,1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	return NO_ERRORS;
}

ERROR BMP_ReadColourPalette(BMPFILE* bmp_fp)
{
	return NO_ERRORS;
}

ERROR BMP_ReadImageData(BMPFILE* bmp_fp)
{
	bmp_fp->imageData = (uint8_t*)malloc(sizeof(uint8_t)*(bmp_fp->bmpInfoHeader->imageSize));
	fseek(bmp_fp->bmpfile_fp,bmp_fp->bmpFileHeader->offset,SEEK_SET);
	if (!(fread((void*)(bmp_fp->imageData),sizeof(uint8_t)*(bmp_fp->bmpInfoHeader->imageSize),1,bmp_fp->bmpfile_fp)==1))
	{
		return FILE_READ_ERROR;
	}
	
	return NO_ERRORS;
}


BMPImage ReadBMP(char* fileName)
{
	BMPImage bmpImage;
	unsigned long int i,j,row,col;
	unsigned long int width, height,size;
	unsigned char r,g,b;
	BMPFILE* bmpfile = CreateBMPFILE(fileName);
	BMP_OpenBitMap(bmpfile,"rb");
	BMP_ReadHeaders(bmpfile);
	
	BMP_ReadImageData(bmpfile);
	
	bmpImage.RGB = (unsigned char**)malloc((bmpfile->bmpInfoHeader->height)*sizeof(unsigned char*));
	for (i=0;i<bmpfile->bmpInfoHeader->height;i++)
	{
		bmpImage.RGB[i] = (unsigned char*)malloc((bmpfile->bmpInfoHeader->width)*sizeof(unsigned char)*3);
	}
	
	width = (unsigned long int)(bmpfile->bmpInfoHeader->width);
	height = (unsigned long int)(bmpfile->bmpInfoHeader->height);
	size = (unsigned long int)(bmpfile->bmpInfoHeader->imageSize);
	
	bmpImage.width = width;
	bmpImage.height = height;
	
	/*copy data*/
	for (row=0;row<height;row++)
	{
		i = row*width*3;
		for (col=0;col<width*3;col+=3)
		{
			j = col;
			b = (unsigned char)(bmpfile->imageData[i+j]);
			g = (unsigned char)(bmpfile->imageData[i+j+1]);
			r = (unsigned char)(bmpfile->imageData[i+j+2]);
			bmpImage.RGB[row][col] = r;
			bmpImage.RGB[row][col+1] = g;
			bmpImage.RGB[row][col+2] = b;
		}
	}
	
	
	BMP_CloseBitMap(bmpfile);
	DestroyBMPFILE(bmpfile);
	
	return bmpImage;
};


