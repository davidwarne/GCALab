/* File: GCALab.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 18/04/2012
 * Last Modified: 04/10/2012
 *
 */

#ifndef __GCALAB_H
#define __GCALAB_H

/*standard headers*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*custom headers*/ 
#include "mesh.h"
#include "GCA.h"
#include "GCALab_fio.h"

#define GCALAB_VERSION 0.4

/*this error code should be consistent with the error codes in mesh.h*/
#define GCALAB_SUCCESS 127
#define GCALAB_MEM_ERROR OUT_OF_MEMORY
#define GCALAB_FATAL_ERROR -128
#define GCALAB_INVALID_OPTION -125
#define GCALAB_UNKNOWN_OPTION -124
#define GCALAB_CL_PARSE_ERROR -122
#define GCALAB_INVALID_WS_ERROR -121;
/*a few new error codes*/

#define GCALAB_GRAPHICS_MODE 0
#define GCALAB_TEXT_MODE 1
#define GCALAB_BATCH_MODE 2

#ifndef GCALAB_DEFAULT_MODE
#define GCALAB_DEFAULT_MODE BATCH_MODE
#endif

#ifndef GCALAB_MAX_WORKSPACES
#define GCALAB_MAX_WORKSPACES 10
#endif

#define GCALAB_MAX_INPUTLENGTH 10

typedef struct CL_Options_struct CL_Options;
typedef struct GCALabworkspace_struct GCALab_WS;


struct GCALabworkspace_struct
{
	GraphCellularAutomaton *GCAList;
	mesh *GCAGeometry;
	unsigned int *commandqueue;
	unsigned int *commandtarget;
	void **commandparams;
	unsigned int numcommands;
	unsigned int cur_command;
	unsigned int state;
}

struct CL_Options_struct
{
	unsigned char mode;
	unsigned char debug;
	/*flags for input and output of CA data*/
	unsigned char load_CA;
	unsigned char save_CA;
	
	char *CAInputFilename;
	char *CAOutputFilename;
	
	/*if not loading then we need to be generating the graph and the CA*/
	unsigned char GenTopology;
	unsigned char GenCA;
	/*topology data if generating*/
	unsigned char type; /*sphere, torus, etc...*/
	unsigned int N;
	unsigned int k;
	
	/*CA definition if generating*/
	unsigned char ECA;
	unsigned char GOL;
	unsigned char numStates;
	unsigned char ruletype;
	unsigned int rule;
	unsigned int windowSize;
	unsigned char ICtype;
	
	/*Simulation info*/
	unsigned int start_t;
	unsigned int end_t;
	unsigned char visualise;
	
	/*analysis to compute*/
	unsigned char do_ShannonE;
	unsigned char do_WordE;
	unsigned char do_InputE;
	unsigned char do_Lambda;
	unsigned char do_Z;
	unsigned char do_G;
	unsigned char do_Pr;
	unsigned char do_R;
	chunk range[2];
	unsigned int numSamples;
};

/*function prototypes*/
char GCALab_Init(int argc,char **argc,CL_Options **opts);
char GCALab_NewWorkSpace(int GCALimit);
char GCALab_QueueCommand(unsigned char ws_id,unsigned int command_id,unsigned int target_id,void *params);
char GCALab_DequeueCommand(unsigned char ws_id,unsigned int index);
char GCALab_ProcessCommandQueue(unsigned char ws_id);
char GCALab_PauseCommandQueue(unsigned char ws_id);
void GCALab_ShutDown(char rc); 
unsigned int GCALab_GetCommandCode(char *cmd);

void GCALab_GraphicsMode(CL_Options* opts);
void GCALab_TextMode(CL_Options* opts);
char **GCALab_CommandPrompt(void);
char GCALab_TestPointer(void* ptr);
char GCALab_PrintWorkSpace(unsigned char ws_id);
char GCALab_ListWorkSpaces(void);
char GCALab_ValidWSId(unsigned int ws_id);
char GCALab_BatchMode(CL_Options* opts);
void GCALab_SplashScreen(void);
void GCALab_HandleErr(char rc);
void PrintGNUGPL3(void);
void PrintAbout(void);
void PrintUsage(void);
#endif
