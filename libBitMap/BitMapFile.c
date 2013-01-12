

#include"BitMapFile.h"

BMPFILE* CreateBMPFILE(char* fileName)
{
	BMPFILE* bmpfile = (BMPFILE*)malloc(sizeof(BMPFILE));
	if (!bmpfile)
	{
		return NULL;
	}
	else
	{
		char c;
		int i,j;
		c = 'a';
		i = 0;
		while (c != '\0')
		{ 
			c = fileName[i];
			i++;
		}
		bmpfile->fileName = (char*)malloc(i);
		for (j=0;j<i;j++)
			bmpfile->fileName[j] = fileName[j]; 
		return bmpfile;
	}
}

BMPFILE* CreateBMPFILE_FromImage(char* fileName,BMPImage* image)
{
	BMPFILE* bmpfile = (BMPFILE*)malloc(sizeof(BMPFILE));
	if (!bmpfile)
	{
		return NULL;	
	}
	else
	{

		char c;
		unsigned long int i,j,k;
		unsigned long int fhSize, ihSize;
		fhSize = 14;
		ihSize = 40;
		c = 'a';
		i = 0;
		while (c != '\0')
		{ 
			c = fileName[i];
			i++;
		}
		bmpfile->fileName = (char*)malloc(i);
		for (j=0;j<i;j++)
			bmpfile->fileName[j] = fileName[j];
		
		/*Allocate memory for headers*/
		
		bmpfile->bmpFileHeader = (BMPFILEHEADER*)malloc(sizeof(BMPFILEHEADER));
		bmpfile->bmpInfoHeader = (BMPINFOHEADER*)malloc(sizeof(BMPINFOHEADER));

		/*now fille the headers*/
		bmpfile->bmpFileHeader->type[0] = 'B';
		bmpfile->bmpFileHeader->type[1] = 'M';
		
		bmpfile->bmpFileHeader->fileSize = (image->width)*(image->height)*3 + fhSize + ihSize;
		bmpfile->bmpFileHeader->reserved1 = 0x1337;
		bmpfile->bmpFileHeader->reserved2 = 0xC0DE;
		bmpfile->bmpFileHeader->offset = fhSize + ihSize;
		
		bmpfile->bmpInfoHeader->size = ihSize;
		bmpfile->bmpInfoHeader->width = (uint32_t)(image->width);
		bmpfile->bmpInfoHeader->height = (uint32_t)(image->height);
		bmpfile->bmpInfoHeader->planes = (uint16_t)1;
		bmpfile->bmpInfoHeader->colorDepth = (uint16_t)24;
		bmpfile->bmpInfoHeader->compression = (uint32_t)0;
		bmpfile->bmpInfoHeader->imageSize = (image->width)*(image->height)*3;
		bmpfile->bmpInfoHeader->hRes = (uint32_t)2835;
		bmpfile->bmpInfoHeader->wRes = (uint32_t)2835;
		bmpfile->bmpInfoHeader->paletteSize = (uint32_t)0;
		bmpfile->bmpInfoHeader->numImportantColours = (uint32_t)0;

		/*allocate memory for BMPFile*/
		bmpfile->imageData = (uint8_t*)malloc(bmpfile->bmpInfoHeader->imageSize);
	
		/*copy data to write*/
		for (i=0;i<image->height;i++)
		{
			k=i*(image->width)*3;
			for (j=0;j<(image->width)*3;j+=3)
			{
				bmpfile->imageData[k+j] = (uint8_t)image->RGB[i][j+2];
				bmpfile->imageData[k+j+1] = (uint8_t)image->RGB[i][j+1];
				bmpfile->imageData[k+j+2] = (uint8_t)image->RGB[i][j];
			}
		}
		
		return bmpfile;	
	}
}

ERROR BMP_OpenBitMap(BMPFILE* bmp_fp, char* option)
{
	/*Open bmp file handle*/
	bmp_fp->bmpfile_fp = fopen(bmp_fp->fileName,option);
	if (!(bmp_fp->bmpfile_fp))
	{
		return FILE_OPEN_ERROR;
	}
	else
	{
		return NO_ERRORS;
	}
}

ERROR BMP_CloseBitMap(BMPFILE* bmp_fp)
{
	if(fclose(bmp_fp->bmpfile_fp))
	{
		return FILE_CLOSE_ERROR;
	}
	
	return NO_ERRORS;
}

ERROR DestroyBMPFILE(BMPFILE* bmp_fp)
{
	free(bmp_fp->fileName);
	free(bmp_fp->bmpFileHeader);
	free(bmp_fp->bmpInfoHeader);
	free(bmp_fp->imageData);
	free(bmp_fp);
	return NO_ERRORS;
}
