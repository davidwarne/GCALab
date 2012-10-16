/* GCALab: An analysis tool for Graph Cellular Automata
 * Copyright (C) 2012  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/* File: GCALab.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 18/04/2012
 * Last Modified: 16/10/2012
 *
 */

#ifndef __GCALAB_H
#define __GCALAB_H

/*standard headers*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*for threading*/
#include <pthread.h>

/*custom headers*/ 
#include "mesh.h"
#include "GCA.h"
#include "GCALab_fio.h"

#define GCALAB_VERSION 0.05

/*this error code should be consistent with the error codes in mesh.h*/
#define GCALAB_SUCCESS 			127
#define GCALAB_MEM_ERROR OUT_OF_MEMORY
#define GCALAB_FATAL_ERROR 		-128
#define GCALAB_INVALID_OPTION 	-125
#define GCALAB_UNKNOWN_OPTION 	-124
#define GCALAB_CL_PARSE_ERROR 	-122
#define GCALAB_INVALID_WS_ERROR -121
#define GCALAB_THREAD_ERROR 	-120
#define GCALAB_QUEUE_FULL	 	-119
/*a few new error codes*/

#define GCALAB_GRAPHICS_MODE 	0
#define GCALAB_TEXT_MODE 		1
#define GCALAB_BATCH_MODE 		2

#ifndef GCALAB_DEFAULT_MODE
#define GCALAB_DEFAULT_MODE GCALAB_TEXT_MODE
#endif

#ifndef GCALAB_MAX_WORKSPACES
#define GCALAB_MAX_WORKSPACES 	10
#endif

#ifndef GCALAB_DEFAULT_MAX_GCA 
#define GCALAB_DEFAULT_MAX_GCA 	25
#endif

#ifndef GCALAB_COMMAND_BUFFER_SIZE
#define GCALAB_COMMAND_BUFFER_SIZE	1024
#endif

#define GCALAB_MAX_STRLEN 128

#ifndef GCALAB_MAX_RESULTS 
#define GCALAB_MAX_RESULTS 2048
#endif


#define GCALAB_WS_STATE_IDLE 		0
#define GCALAB_WS_STATE_PROCESSING 	1
#define GCALAB_WS_STATE_PAUSED 		2
#define GCALAB_WS_STATE_EXITING 	3
#define GCALAB_WS_STATE_ERROR 		-1

#define GCALAB_NOP 		0
#define GCALAB_LOAD 	1
#define GCALAB_SAVE 	2
#define GCALAB_SIMULATE 3
#define GCALAB_GCA 		4
#define GCALAB_ENTROPY 	5
#define GCALAB_PARAM	6
#define GCALAB_REVERSE	7

#define WS(a) GCALab_Global[(a)]

typedef struct GCALab_CL_Options_struct GCALab_CL_Options;
typedef struct GCALabworkspace_struct GCALab_WS;
typedef struct GCALabOutput_struct GCALabOutput;


struct GCALabOutput_struct
{
	int type;
	char id[GCALAB_MAX_STRLEN];
	unsigned int datalen;
	void *data;
};

struct GCALabworkspace_struct
{
	GraphCellularAutomaton **GCAList;
	mesh **GCAGeometry;
	int numGCA;
	int maxGCA;
	int numresults;
	GCALabOutput **results;
	unsigned int *	commandqueue;
	unsigned int *commandtarget;
	char ***commandparams;
	int *numparams;
	unsigned int numcommands;
	unsigned int qhead;
	unsigned int qtail;
	unsigned int state;
	pthread_t worker;
	pthread_mutex_t wslock;
};

struct GCALab_CL_Options_struct
{
	unsigned char mode;
	unsigned char debug;
	/*flags for input and output of CA data*/
	unsigned char load_CA;
	unsigned char save_CA;
	
	char *CAInputFilename;
	char *CAOutputFilename;
};

/*function prototypes*/
char GCALab_Init(int argc,char **argv,GCALab_CL_Options **opts);
char GCALab_NewWorkSpace(int GCALimit);
void GCALab_LockWS(unsigned char ws_id);
void GCALab_UnLockWS(unsigned char ws_id);
char GCALab_QueueCommand(unsigned char ws_id,unsigned int command_id,unsigned int target_id,char **params,int numparams);
char GCALab_CancelCommand(unsigned char ws_id,unsigned int index);
char GCALab_ProcessCommandQueue(unsigned char ws_id);
char GCALab_PauseCommandQueue(unsigned char ws_id);
char GCALab_PrintCommandQueue(unsigned char ws_id);

unsigned int GCALab_GetState(unsigned char ws_id);
void GCALab_SetState(unsigned char ws_id,unsigned int state);
void GCALab_ShutDown(char rc); 
unsigned int GCALab_GetCommandCode(char *cmd);

void *GCALab_Worker(void *params);
char GCALab_DoNextCommand(unsigned char ws_id);
unsigned char GCALab_CommandQueueEmpty(unsigned char ws_id);

void GCALab_GraphicsMode(GCALab_CL_Options* opts);
void GCALab_TextMode(GCALab_CL_Options* opts);
char GCALab_BatchMode(GCALab_CL_Options* opts);

char **GCALab_CommandPrompt(int *numargs);
char GCALab_PrintWorkSpace(unsigned char ws_id);
char GCALab_ListWorkSpaces(void);
void GCALab_PrintHelp(void);

void GCALab_SplashScreen(void);
void GCALab_PrintLicense(void);
void GCALab_PrintAbout(void);
void GCALab_PrintUsage(void);

void GCALab_HandleErr(char rc);
char GCALab_TestPointer(void* ptr);
char GCALab_ValidWSId(unsigned int ws_id);

void GCALab_InitCL_Options(GCALab_CL_Options* opts);
GCALab_CL_Options* GCALab_ParseCommandLineArgs(int argc, char **argv);
char** strvncpy(char **strv,int c,int n);
#endif
