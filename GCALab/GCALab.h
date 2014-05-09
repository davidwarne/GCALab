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
 * Last Modified: 18/01/2013
 *
 */
/**
 * @file GCALab.h
 *
 * @brief
 * @details
 * 
 * @author
 *
 * @version
 * @date
 * @copyright
 */

#ifndef __GCALAB_H
#define __GCALAB_H

/*standard headers*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
/*for threading*/
#include <pthread.h>

#ifdef WITH_GRAPHICS
#include <GL/glut.h>	
#define GCALAB_GRAPHICS_INTERACTION_NONE 0
#define GCALAB_GRAPHICS_INTERACTION_ROTATE 1
#define GCALAB_GRAPHICS_INTERACTION_TRANSLATE 2
#define GCALAB_GRAPHICS_INTERACTION_ZOOM 3
#define GCALAB_GRAPHICS_OFFSET 2
#endif

/*custom headers*/ 
#include "mesh.h"
#include "GCA.h"
#include "GCALab_fio.h"


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
#define GCALAB_MAX_WORKSPACES 	20
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
#define GCALAB_WS_STATE_ERROR 		5

#define GCALAB_NOP 		0
#define GCALAB_LOAD 	1
#define GCALAB_SAVE 	2
#define GCALAB_SIMULATE 3
#define GCALAB_GCA 		4
#define GCALAB_ENTROPY 	5
#define GCALAB_PARAM	6
#define GCALAB_REVERSE	7
#define GCALAB_STATE_FREQUENCIES 8

#define GCALAB_SHANNON_ENTROPY 	0
#define GCALAB_WORD_ENTROPY 	1
#define GCALAB_INPUT_ENTROPY 	2
#define GCALAB_ALL_ENTROPY 		3

#define GCALAB_LAMBDA_PARAM		0
#define GCALAB_Z_PARAM 			1
#define GCALAB_G_PARAM			2
#define GCALAB_C_PARAM			3
#define GCALAB_T_PARAM			4
#define GCALAB_ALL_PARAM		5

#define WS(a) GCALab_Global[(a)]

#ifndef GCALAB_MAXNUM_CMDS
#define GCALAB_MAXNUM_CMDS 14
#endif

#ifndef GCALAB_MAXNUM_OPS
#define GCALAB_MAXNUM_OPS 20
#endif

typedef struct GCALab_CL_Options_struct GCALab_CL_Options;
typedef struct GCALabworkspace_struct GCALab_WS;
typedef struct GCALabOutput_struct GCALabOutput;
typedef struct GCALab_Command_struct GCALab_Cmd;
typedef struct GCALab_Operation_struct GCALab_Op;

/*output data from compute operations*/
struct GCALabOutput_struct
{
	int type;
	char id[GCALAB_MAX_STRLEN];
	unsigned int datalen;
	void *data;
};

/*memory and command queue attached to a thread*/
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

/*command line options*/
struct GCALab_CL_Options_struct
{
	unsigned char mode;
	unsigned char debug;
	/*flags for input and output of CA data*/
	unsigned char load_CA;
	unsigned char save_CA;
	char *ScriptFile;	
	char *CAInputFilename;
	char *CAOutputFilename;
};

/*high level GCALab commands - workspace level*/
struct GCALab_Command_struct
{
	/*the command name*/
	char *id;
	/*function pointer*/
	char (*f)(int,char**);
	/*help infomation*/
	char *args;
	char *desc;
};

/*low level compute operations - CA level*/
struct GCALab_Operation_struct
{
	/*the operation name*/
	char *id;
	/*function pointer*/
	char (*f)(unsigned char, unsigned int,int,char**,GCALabOutput **res);
	/*help infomation*/
	char *args;
	char *desc;
};

/*function prototypes*/
char GCALab_Init(int argc,char **argv,GCALab_CL_Options **opts);
void GCALab_Register_Command(char *id,char (*f)(int, char**),char * args, char * desc);
void GCALab_Register_Operation(char *id,char (*f)(unsigned char,unsigned int,int,char**,GCALabOutput**),char* args,char * desc);
char GCALab_Process_Command(int nargs,char ** args);
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

/*Graphics Mode*/

#ifdef WITH_GRAPHICS
void GCALab_Graphics_Init(void);
void GCALab_Graphics_Display(void);
void GCALab_Graphics_KeyPressed (unsigned char key, int x, int y);
void GCALab_Graphics_SpecialKeyPressed(int key, int x, int y);
void GCALab_Graphics_MouseClick(int button,int state,int x, int y);
void GCALab_Graphics_MouseMove(int x, int y);
void GCALab_Graphics_Reshape(int w, int h);
void GCALab_Graphics_Timer(int x);
void GCALab_Graphics_DrawGCA(unsigned int ws_id,unsigned int gca_id);
void GCALab_Graphics_DrawBBox(void);
void GCALab_Graphics_DrawGrid(void);
#endif
/*Text and Batch Mode*/
char **GCALab_CommandPrompt(int *numargs);
char GCALab_PrintWorkSpace(unsigned char ws_id);
char GCALab_ListWorkSpaces(void);
void GCALab_PrintHelp(void);
void GCALab_PrintOperations(void);
char GCALab_PrintCA(unsigned char ws_id,unsigned int gca_id);
char GCALab_PrintSTP(unsigned char ws_id,unsigned int gca_id);
char GCALab_PrintResults(unsigned char ws_id,unsigned int res_id);

char ** GCALab_ReadScriptCommand(char * filename, int *numargs);

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


/*menu commands*/
char GCALab_CMD_NewWorkSpace(int argc, char **argv);
char GCALab_CMD_PrintWorkSpace(int argc, char **argv);
char GCALab_CMD_ListWorkSpaces(int argc, char **argv);
char GCALab_CMD_PrintHelp(int argc, char **argv);
char GCALab_CMD_PrintOperations(int argc,char **argv);
char GCALab_CMD_ChangeWorkSpace(int argc, char **argv);
char GCALab_CMD_QueueCommand(int argc, char **argv);
char GCALab_CMD_DeleteCommand(int argc, char **argv);
char GCALab_CMD_ExecuteQueue(int argc, char **argv);
char GCALab_CMD_StopQueue(int argc, char **argv);
char GCALab_CMD_PrintCA(int argc, char **argv);
char GCALab_CMD_PrintSTP(int argc, char ** argv);
char GCALab_CMD_PrintResults(int argc, char **argv);
char GCALab_CMD_Quit(int argc, char **argv);

/*compute operations*/
char GCALab_OP_NOP(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Load(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Save(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Simulate(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_GCA(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Entropy(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Param(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Reverse(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Freq(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
char GCALab_OP_Pop(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res);
#endif
#endif
