/* File: GCALab_fio.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 19/04/2012
 * Last Modified: 19/04/2012
 *
 * Description: File IO function defintions for GCALab formats
 *
 *==============================================================================
 */
 
#ifndef __GCALAB_FIO_H
#define __GCALAB_FIO_H

#include <stdio.h>
#include "GCA.h"

#define FLOAT32 0
#define FLOAT64 2
#define UINT8 3
#define UINT16 4
#define UINT32 5
#define UINT64 6
#define SINT8 7
#define SINT16 8
#define SINT32 9
#define SINT64 10 
#define HEX8 11
#define HEX16 12
#define HEX32 13
#define HEX64 14

#if CHUNK_SIZE_BITS == 64
	#define CHUNK HEX64
#elif CHUNK_SIZE_BITS == 32
	#define CHUNK HEX32
#elif CHUNK_SIZE_BITS == 16
	#define CHUNK HEX16
#elif CHUNK_SIZE_BITS == 8
	#define CHUNK HEX8
#endif

/*for now this is it*/
char GCALab_fio_saveCA(char *filename,GraphCellularAutomaton *GCA);

char GCALab_fio_appendData(char* filename,char * name, void * data, int N,unsigned char type);

#endif
