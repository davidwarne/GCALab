/* File: BitMapWriter.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 19/04/2010
 *
 * Summary: Implementations of BitMap Writer Functions.
 *
 */
 
#include "BitMapWriter.h"

ERROR BMP_WriteHeaders(BMPFILE* bmpfile)
{
	if (!(fwrite((void*)(bmpfile->bmpFileHeader->type),2,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpFileHeader->fileSize),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpFileHeader->reserved1),2,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpFileHeader->reserved2),2,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpFileHeader->offset),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->size),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->width),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->height),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->planes),2,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->colorDepth),2,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->compression),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->imageSize),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->hRes),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->wRes),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->paletteSize),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	if (!(fwrite((void*)&(bmpfile->bmpInfoHeader->numImportantColours),4,1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	return NO_ERRORS;
}

ERROR BMP_WriteColourPalette(BMPFILE* bmpfile)
{
	return NO_ERRORS;
}

ERROR BMP_WriteImageData(BMPFILE* bmpfile)
{
	if(!(fwrite((void*)(bmpfile->imageData),(bmpfile->bmpInfoHeader->imageSize),1,bmpfile->bmpfile_fp)==1))
	{
		return FILE_WRITE_ERROR;
	}
	return NO_ERRORS;
}
 


ERROR WriteBMP(char* fileName,BMPImage* image)
{
	ERROR error;
	BMPFILE* bmpfile = CreateBMPFILE_FromImage(fileName,image);
	error = BMP_OpenBitMap(bmpfile,"wb");
	error = BMP_WriteHeaders(bmpfile);
	error = BMP_WriteImageData(bmpfile);
	error = BMP_CloseBitMap(bmpfile);
	error = DestroyBMPFILE(bmpfile);
	return error;
}
