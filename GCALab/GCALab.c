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
/* File: GCALab.c
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
 * Version History:
 *       v 0.01 (18/04/2012) - i. Initial Version...  was bored at work and 
 *                             thought I should write a main program for my research.
 *       v 0.02 (19/04/2012) - i. Implemented and tested Commandline parsing
 *                             ii. starting to implement Batch mode
 *       v 0.04 (04/10/2012) - i. Re-implementing main processing loop entirely
 *                                this is a major restructure.
 *       v 0.05 (15/10/2012) - i. implemented the command prompt,
 *                             ii. implemented thread creation and deletion
 *                             iii. implemented new workspace command
 *       v 0.06 (16/10/2012) - i. commands for queue management implemented.
 *                             ii. Began work on main compute routine.
 *                             iii. Removed legacy code
 *       v 0.07 (18/10/2012) - i. Implemented gca creation commands
 *                             ii. Implemented Entropy commands
 *                             iii. Implemented parameter commands
 *                             iv. Implement pre-image comput command
 *       v 0.08 (19/10/2012) - i. Testing of compute functions
 *                             ii. save function implemented
 *       v 0.09 (20/10/2012) - i. load function implemented
 *                             ii. simulation implemented 
 *                             iii. testing of functions for no-eca gca
 *       v 0.10 (21/10/2012) - i. implemented batch mode
 *                             ii. implemented State frequency calculation
 *       v 0.11 (03/11/2012) - i. created command and operation structs
 *                                which use function pointers to made the software easiliy
 *                                extendable.
 *                             ii. re-impement text and batch command using this new
 *                                 function pointer framework.
 *                             iii. Impemented compute operations using this new function
 *                                  pointer framework.
 *                             iv. removed switch statement from DoNextCommand and replaced
 *                                  it with a call to a function pointer.
 *       v 0.12 (14/12/2012) - i. Implemented Initial version OpenGL callbacks
 *                             ii. basic CA animation now supported in Display func.
 *                             iii. made text mode commands accessable via 'c' key
 *       v 0.13 (11/01/2013) - i. Happy new year!
 *                             ii. Added mouse handlers
 *                             iii. cleaned up the display function.
 *       v 0.14 (12/01/2013) - i. Fixed but in parsing rule code, should be and unsigned int
 *                                not an unsigned char.
 *                             ii. added synced version of CA evolution, good for animations.
 *       v 0.15 (13/01/2013) - i. improved the graphics quality
 *                             ii. added bitmap output for time-space patterns.
 *       v 0.16 (16/01/2013) - i. added lut coloring option in graphics mode
 *                             ii. fixed a bug in gca create command when creating 
 *                                 an ECA with a user specified window size.
 *                             iii. added a all parameters option for the param
 *                                  command.
 *
 * Description: Main Program for Graph Cellular Automata generation, simulation,
 *              analysis and Visualisation.
 *
 * TODO List:
 *		1. test commandline parsing - done (v 0.02)
 *		2. implement batch mode first - cancelled (v 0.04)
 *		3. use ptheads to implement a main loop plus compute threads - done (v 0.05)
 * 		4. implement text mode first, then re-implement batch mode, then graphics 
 * 		5. carefully separated main engine from the mode type
 * 		6. remember to comment function pointer framework carefully.
 * 		7. should add a cycle compute function, possibly useful comparison
 * 		8. build in more plotting/image saving options
 * 		9. add a create animation option
 * Known Issues:
 *
 *==============================================================================
 */
#define GCALAB_VERSION 0.14
#include "GCALab.h"
#ifdef __linux__
/*global array of workspace addresses*/
GCALab_WS **GCALab_Global;
unsigned int GCALab_numWS;
GCALab_Cmd GCALab_Cmds[14];
unsigned int GCALab_numCmds;
GCALab_Op GCALab_Ops[9];
unsigned int GCALab_numOps;
unsigned char GCALab_mode;
unsigned int cur_ws;
#ifdef WITH_GRAPHICS
/* light settings*/
GLfloat ambientLight[] = { 0.1f, 0.1f, 0.1f, 1.0f };
GLfloat diffuseLight[] = { 3.0f, 3.0f, 3.0f, 1.0f };
GLfloat position[] = { 2.0f, 2.0f, 2.0f, 1.0f };
GLfloat positionMirror[] = { 2.0f, -2.0f, 2.0f, 1.0f };
GLfloat fogCol[] = {0.5f,0.5f,0.5f,1.0f};
GLUquadric* quad;
float translateX,translateY,zoom;
float theta, phi;
unsigned int cur_gca;
unsigned int cur_res;
unsigned char moving;
int xprev;
int yprev;
float cellColf[8][3] = {{0,0,0},{1,1,1},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1}};
float lutColf[16][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1},{1,1,1},
						{1,0.5,0},{0.0,1,0.5},{0.5,0,1},{1,0.5,1},{1,1,0.5},{0.5,1,1},{0.5,0.5,0.5}};
unsigned char lutColmode;
#endif
char stateInitials[6] = {'I','R','P','Q','E'};
char* statenames[6] = {"Idle","Running","Paused","Exiting","Error"};
char cellsymbols[8] = {' ','O','*','-','x','+','#','^'};
/* main(): entry point of the GCALab Program
 */
int main(int argc, char **argv)
{
	GCALab_CL_Options* CL_opt;
	char rc;
#ifdef WITH_GRAPHICS
	glutInit(&argc,argv);
#endif
	rc = GCALab_Init(argc,argv,&CL_opt);
	GCALab_HandleErr(rc);
	/*Start up in either batch or interactive mode*/
	switch(GCALab_mode)
	{
		case GCALAB_GRAPHICS_MODE:
			GCALab_GraphicsMode(CL_opt);
			break;
		case GCALAB_TEXT_MODE:
			GCALab_TextMode(CL_opt);
			break;
		case GCALAB_BATCH_MODE:
			GCALab_BatchMode(CL_opt);
			break;
		default:
			fprintf(stderr,"ERROR: You Should not see this!\n");
			exit(1);
			break;
	}
	
	free(CL_opt);
	pthread_exit(NULL);
} 

/* GCALab_TextMode(): starts a GCALab session in text-only mode
 */
void GCALab_TextMode(GCALab_CL_Options* opts)
{
	char **userinput;
	int numargs;
	char* cmd;
	char rc;
	GCALab_SplashScreen();
	
	while(1)
	{
		int i;
		/*get user input*/
		userinput = GCALab_CommandPrompt(&numargs);
		/*ensure input is valid*/
		rc = GCALab_TestPointer((void*)userinput);
		GCALab_HandleErr(rc);
		
		cmd = userinput[0];
		for (i=0;i<GCALab_numCmds;i++)
		{
			if (!strcmp(cmd,GCALab_Cmds[i].id))
			{
				rc = (*(GCALab_Cmds[i].f))(numargs,userinput);
				GCALab_HandleErr(rc);
			}
		}
		/*clean up user args*/
		for (i=0;i<numargs;i++)
		{
			free(userinput[i]);
		}
		free(userinput);
	}
	return;
}

/* GCALab_InteractiveMode(): Runs GCALab interactively
 */
void GCALab_GraphicsMode(GCALab_CL_Options* opts)
{
#ifdef WITH_GRAPHICS
	int hPix,wPix;
	GCALab_SplashScreen();
	/*init display mode*/
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	/*get the screen resolution*/
	hPix = glutGet(GLUT_SCREEN_HEIGHT);
	wPix = glutGet(GLUT_SCREEN_WIDTH);
	/*set up window*/
	glutInitWindowSize(hPix,hPix);
	glutInitWindowPosition(wPix - hPix,0);
	glutCreateWindow("GCALab");
	GCALab_Graphics_Init();
	glutDisplayFunc(GCALab_Graphics_Display);
	glutReshapeFunc(GCALab_Graphics_Reshape);
	glutKeyboardFunc(GCALab_Graphics_KeyPressed);
	glutSpecialFunc(GCALab_Graphics_SpecialKeyPressed);
	glutMouseFunc(GCALab_Graphics_MouseClick);
	glutMotionFunc(GCALab_Graphics_MouseMove);
	glutTimerFunc(1000/60,GCALab_Graphics_Timer,1);

	cur_gca = 0;
	cur_res = 0;
	lutColmode = 0;
	moving = GCALAB_GRAPHICS_INTERACTION_NONE;	
	glutMainLoop();
#else
	fprintf(stderr,"Graphics Mode not enabled in this build.\n");
#endif
	return;
}

/* GCALab_BatchMode(): Runs Commands is batch mode, 
 *                     i.e., reads commands from a file
 */
char GCALab_BatchMode(GCALab_CL_Options* opts)
{
	char **userinput;
	int numargs;
	char* cmd;
	char rc;
	GCALab_SplashScreen();
	
	while(1)
	{
		int i;
		/*get user input*/
		userinput = GCALab_ReadScriptCommand(opts->ScriptFile,&numargs);
		opts->ScriptFile = NULL;
		/*ensure input is valid*/
		rc = GCALab_TestPointer((void*)userinput);
		GCALab_HandleErr(rc);
		
		cmd = userinput[0];
		for (i=0;i<GCALab_numCmds;i++)
		{
			if (!strcmp(cmd,GCALab_Cmds[i].id))
			{
				rc = (*(GCALab_Cmds[i].f))(numargs,userinput);
				GCALab_HandleErr(rc);
			}
		}

		/*clean up user args*/
		for (i=0;i<numargs;i++)
		{
			free(userinput[i]);
		}
		free(userinput);
	}
}

/* GCALab_TestPointer(): returns an error if the pointer is null
 */
char GCALab_TestPointer(void* ptr)
{
	if (!ptr)
		return GCALAB_MEM_ERROR;
	else
		return GCALAB_SUCCESS;
}

/* GCALab_ValidWSId(): tests if the given workspace id maps to a real workspace
 */
char GCALab_ValidWSId(unsigned int ws_id)
{
	/*for now simple*/
	if (ws_id >= GCALab_numWS)
		return GCALAB_INVALID_WS_ERROR;
	else
		return GCALAB_SUCCESS;
}

/* GCALab_HandleErr(): tests if the given returen code indicates an error,
 *                     if the error is fatal then it aborts.
 */
void GCALab_HandleErr(char rc)
{
	if (rc <= 0)
	{
		switch (rc)
		{
			case GCALAB_MEM_ERROR:
				fprintf(stderr,"OUT OF MEMORY!!!\n");
				exit(1);
				break;
			case GCALAB_FATAL_ERROR:
				fprintf(stderr,"A fatal error occurred. Aborting.\n");
				exit(1);
				break;
			case GCALAB_INVALID_OPTION:
				break;
			case GCALAB_UNKNOWN_OPTION:
				break;
			case GCALAB_CL_PARSE_ERROR:
				break;
			case GCALAB_INVALID_WS_ERROR:
				break;
			case GCALAB_THREAD_ERROR:
				exit(1);
				break;
		}
	}
}

/* GCALab_Init(): sets up GCALab with default settings, unless
 *                overridden via start-up commands
 */
char GCALab_Init(int argc,char **argv,GCALab_CL_Options **opts)
{
	GCALab_Global = (GCALab_WS **)malloc(GCALAB_MAX_WORKSPACES*sizeof(GCALab_WS*));
	if(!(GCALab_Global))
	{
		return GCALAB_FATAL_ERROR;
	}

	/*register commands*/
	GCALab_numCmds = 14;
	GCALab_Cmds[0].id = "new-work";
	GCALab_Cmds[0].f = &GCALab_CMD_NewWorkSpace;
	GCALab_Cmds[0].args = "n";
	GCALab_Cmds[0].desc = "Creates a new workspace with n CA slots.";

	GCALab_Cmds[1].id = "print-work";
	GCALab_Cmds[1].f = &GCALab_CMD_PrintWorkSpace;
	GCALab_Cmds[1].args = "none";
	GCALab_Cmds[1].desc = "Print the current workspace.";

	GCALab_Cmds[2].id = "list-all-work";
	GCALab_Cmds[2].f = &GCALab_CMD_ListWorkSpaces;
	GCALab_Cmds[2].args = "none";
	GCALab_Cmds[2].desc = "Print Summary of all workspaces.";

	GCALab_Cmds[3].id = "ch-work";
	GCALab_Cmds[3].f = &GCALab_CMD_ChangeWorkSpace;
	GCALab_Cmds[3].args = "id";
	GCALab_Cmds[3].desc = "Changes current workspace to be id.";
	
	GCALab_Cmds[4].id = "q-cmd";
	GCALab_Cmds[4].f = &GCALab_CMD_QueueCommand;
	GCALab_Cmds[4].args = "cmd";
	GCALab_Cmds[4].desc = "Enqueues the GCA operation to the current command queue.";
	
	GCALab_Cmds[5].id = "del-cmd";
	GCALab_Cmds[5].f = &GCALab_CMD_DeleteCommand;
	GCALab_Cmds[5].args = "cmd_id";
	GCALab_Cmds[5].desc = "Sets GCA operation to cmd_id to be ignored.";
	
	GCALab_Cmds[6].id = "exec-q";
	GCALab_Cmds[6].f = &GCALab_CMD_ExecuteQueue;
	GCALab_Cmds[6].args = "none";
	GCALab_Cmds[6].desc = "The current queue will start processing.";
	
	GCALab_Cmds[7].id = "stop-q";
	GCALab_Cmds[7].f = &GCALab_CMD_StopQueue;
	GCALab_Cmds[7].args = "none";
	GCALab_Cmds[7].desc = "The current queue will pause after completion of the current operation.";
	
	GCALab_Cmds[8].id = "print-ca";
	GCALab_Cmds[8].f = &GCALab_CMD_PrintCA;
	GCALab_Cmds[8].args = "ca_id";
	GCALab_Cmds[8].desc = "Prints GCA with id ca_id in the current workspace.";
	
	GCALab_Cmds[9].id = "print-st";
	GCALab_Cmds[9].f = &GCALab_CMD_PrintSTP;
	GCALab_Cmds[9].args = "ca_id";
	GCALab_Cmds[9].desc = "Prints GCA evolution.";
	
	GCALab_Cmds[10].id = "print-result";
	GCALab_Cmds[10].f = &GCALab_CMD_PrintResults;
	GCALab_Cmds[10].args = "res_id";
	GCALab_Cmds[10].desc = "Prints result data with id res_id.";
	
	GCALab_Cmds[11].id = "quit";
	GCALab_Cmds[11].f = &GCALab_CMD_Quit;
	GCALab_Cmds[11].args = "none";
	GCALab_Cmds[11].desc = "Exit GCALab.";
 	
	GCALab_Cmds[12].id = "help";
	GCALab_Cmds[12].f = &GCALab_CMD_PrintHelp;
	GCALab_Cmds[12].args = "none";
	GCALab_Cmds[12].desc = "Prints this help menu.";
	
	GCALab_Cmds[13].id = "list-cmds";
	GCALab_Cmds[13].f = &GCALab_CMD_PrintOperations;
	GCALab_Cmds[13].args = "none";
	GCALab_Cmds[13].desc = "Prints available compute commands.";
	/*register operations*/	
	GCALab_Ops[GCALAB_NOP].id = "nop";
	GCALab_Ops[GCALAB_NOP].f = &GCALab_OP_NOP;
	GCALab_Ops[GCALAB_NOP].args = "none";
	GCALab_Ops[GCALAB_NOP].desc = "No Operation";
	GCALab_Ops[GCALAB_LOAD].id = "load";
	GCALab_Ops[GCALAB_LOAD].f = &GCALab_OP_Load;
	GCALab_Ops[GCALAB_LOAD].args = "i -f filename";
	GCALab_Ops[GCALAB_LOAD].desc = "Loads a *.gca file into the current workspace";
	GCALab_Ops[GCALAB_SAVE].id = "save";
	GCALab_Ops[GCALAB_SAVE].f = &GCALab_OP_Save;
	GCALab_Ops[GCALAB_SAVE].args = "i -f filename (-g | -r)";
	GCALab_Ops[GCALAB_SAVE].desc = "Save a GCA or result with id to file";
	GCALab_Ops[GCALAB_SIMULATE].id = "simulate";
	GCALab_Ops[GCALAB_SIMULATE].f = &GCALab_OP_Simulate;
	GCALab_Ops[GCALAB_SIMULATE].args = "i -t Tfinal [-I] [-f icfile | -c ictype]";
	GCALab_Ops[GCALAB_SIMULATE].desc = "simulates the id to Tfinal";
	GCALab_Ops[GCALAB_GCA].id = "gca";
	GCALab_Ops[GCALAB_GCA].f = &GCALab_OP_GCA;
	GCALab_Ops[GCALAB_GCA].args = "i (((-m meshfile | -t numcells genus) -s numstates -r ruletype rulecode) | -eca numCells numNeighbours rule) [-c ictype] [-w windowsize]";
	GCALab_Ops[GCALAB_GCA].desc = "Creates a new graph cellular automaton in the current workspace";
	GCALab_Ops[GCALAB_ENTROPY].id = "entropy";
	GCALab_Ops[GCALAB_ENTROPY].f = &GCALab_OP_Entropy;
	GCALab_Ops[GCALAB_ENTROPY].args = "i -n numsamples -t timesteps -e entropytype";
	GCALab_Ops[GCALAB_ENTROPY].desc = "Computes entropy measures of graph cellular automaton at i";
	GCALab_Ops[GCALAB_PARAM].id = "param";
	GCALab_Ops[GCALAB_PARAM].f = &GCALab_OP_Param;
	GCALab_Ops[GCALAB_PARAM].args = "i -p paramtype [-l config0 configN]";
	GCALab_Ops[GCALAB_PARAM].desc = "Computes complexity parameters such as Langton's lambda";
	GCALab_Ops[GCALAB_REVERSE].id = "reverse";
	GCALab_Ops[GCALAB_REVERSE].f = &GCALab_OP_Reverse;
	GCALab_Ops[GCALAB_REVERSE].args = "i";
	GCALab_Ops[GCALAB_REVERSE].desc = "Computes pre-images of the current configuration of the graph cellular automaton at i";
	GCALab_Ops[GCALAB_STATE_FREQUENCIES].id = "freq";
	GCALab_Ops[GCALAB_STATE_FREQUENCIES].f = &GCALab_OP_Freq;
	GCALab_Ops[GCALAB_STATE_FREQUENCIES].args = "i (-n numsamples | -l config0 configN)";
	GCALab_Ops[GCALAB_STATE_FREQUENCIES].desc = "Computes state frequency histogram for each cell in the graph cellular automaton at i";
	GCALab_numOps = 9;
	opts[0] = GCALab_ParseCommandLineArgs(argc,argv);
	if (!(opts[0]))
	{
		return GCALAB_CL_PARSE_ERROR;
	}
	GCALab_numWS = 0;
	GCALab_mode = (opts[0])->mode;
	return GCALAB_SUCCESS;
}

/* GCALab_GetCommandCode(): converts the command string into a
 *                          command id
 */
unsigned int GCALab_GetCommandCode(char* cmd)
{
	unsigned int i;

	for (i=0;i<GCALab_numOps;i++)
	{
		if(!strcmp(cmd,GCALab_Ops[i].id)) return i;
	}
	return GCALAB_NOP;
}

/* GCALab_Worker(): Processing thread attached to a workspace
 */
void * GCALab_Worker(void *params)
{
	unsigned char ws_id;
	unsigned char exit,error;
	char rc;
	ws_id = (unsigned int)(unsigned long long)params;

	exit = 0;
	error = 0;
	while(!exit)
	{
		switch(GCALab_GetState(ws_id))
		{
			case GCALAB_WS_STATE_ERROR:
				exit = 1;
				error = 1;
				break;
			case GCALAB_WS_STATE_PAUSED:
				/*just rest for a while*/
				sleep(1);
				break;
			case GCALAB_WS_STATE_PROCESSING:
				/*the next command*/
				rc = GCALab_DoNextCommand(ws_id);
				if (rc <=0 ) GCALab_SetState(ws_id,GCALAB_WS_STATE_ERROR);
				if (GCALab_GetState(ws_id) == GCALAB_WS_STATE_PROCESSING)
				{
					GCALab_SetState(ws_id,GCALAB_WS_STATE_IDLE);
				}
				break;
			case GCALAB_WS_STATE_IDLE:
				/*check if a command is available*/
				if(!GCALab_CommandQueueEmpty(ws_id))
				{
					/*if so then set state to processing*/
					GCALab_SetState(ws_id,GCALAB_WS_STATE_PROCESSING);
				}
				break;
			case GCALAB_WS_STATE_EXITING:
				exit = 1;
				break;
		}
		usleep(1);
	}
	pthread_exit((void*)((long long)error));
}

/* GCALab_CommandQueueEmpty(): evaluates to true when the Command queue for 
 *                             the given workspace is empty.
 */
unsigned char GCALab_CommandQueueEmpty(unsigned char ws_id)
{
	unsigned char empty;
	GCALab_LockWS(ws_id);
	empty = (WS(ws_id)->numcommands == 0);
	GCALab_UnLockWS(ws_id);
	return empty;
}

/* GCALab_DoNextCommand(): executes the next command in the given command queue.
 */
char GCALab_DoNextCommand(unsigned char ws_id)
{
	unsigned int cmd_id,trgt_id;
	int nparams;
	int index,i;
	char **params;
	char rc;
	GCALabOutput *res;

	res = NULL;
	
	GCALab_LockWS(ws_id);
	/*get command data*/
	index = WS(ws_id)->qhead;
	cmd_id = WS(ws_id)->commandqueue[index];
	trgt_id = WS(ws_id)->commandtarget[index];
	nparams = WS(ws_id)->numparams[index];
	params = strvncpy(WS(ws_id)->commandparams[index],nparams,GCALAB_MAX_STRLEN);
	/*remove command from queue*/
	WS(ws_id)->qhead = (WS(ws_id)->qhead+1)%GCALAB_COMMAND_BUFFER_SIZE;
	WS(ws_id)->numcommands--;
	if (nparams != 0)
	{
		WS(ws_id)->numparams[index] = 0;
		for(i=0;i<nparams;i++)
		{
			free(WS(ws_id)->commandparams[index][i]);
		}
		free(WS(ws_id)->commandparams[index]);
	}
	GCALab_UnLockWS(ws_id);

	/*do processing*/
	rc = (*(GCALab_Ops[cmd_id].f))(ws_id,trgt_id,nparams,params, &res);

	/*clean up copied params*/
	for (i=0;i<nparams;i++)
	{
		free(params[i]);
	}
	free(params);
	
	GCALab_LockWS(ws_id);
	if (res != NULL)
	{
		if (WS(ws_id)->numresults < GCALAB_MAX_RESULTS)
		{
			WS(ws_id)->results[WS(ws_id)->numresults] = res;
			WS(ws_id)->numresults++;
		}
	}
	/*Cleanup command*/
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}
/* PrintsAbout(): Prints Author and affiliation information
 */
void GCALab_PrintAbout(void)
{
	fprintf(stdout,"\n The Graph Cellular Automata Lab is a multi-threaded,\n");
	fprintf(stdout," flexible, analysis tool forthe exploration of cellular automata\n"); 
	fprintf(stdout," defined on graphs.\n");
	fprintf(stdout,"\n Version:\t%0.2f\n",GCALAB_VERSION);
	fprintf(stdout," Author:\tDavid J. Warne\n");
	fprintf(stdout," School:\tSchool of Electrical Engineering and Computer Science,\n");
	fprintf(stdout," \t\tThe Queensland University of Technology, Brisbane, Australia\n");
	fprintf(stdout," Contact:\tdavid.warne@qut.edu.au\n");
	fprintf(stdout," Website:\tcoming soon!\n");
	return;
}

/* GCALab_PrintLicense(): prints user license summary
 */
void GCALab_PrintLicense(void)
{
	fprintf(stdout,"\n GCALab v %0.2f  Copyright (C) 2012  David J. Warne\n",GCALAB_VERSION);
	fprintf(stdout," This program comes with ABSOLUTELY NO WARRANTY.\n");
	fprintf(stdout," This is free software, and you are welcome to redistribute it\n");
	fprintf(stdout," under certain conditions.\n");
	return;
}

/* GCALab_SplashScreen(): startup screen
 */
void GCALab_SplashScreen(void)
{
	switch(GCALab_mode)
	{
		case GCALAB_TEXT_MODE:
			fprintf(stdout,"\n +++ Welcome to Graph Cellular Automata Lab! +++\n");
			GCALab_PrintAbout();
			GCALab_PrintLicense();
			fprintf(stdout,"\n Type \'help\' to get started.\n");
			break;
		case GCALAB_BATCH_MODE:
			GCALab_PrintLicense();
			break;
		case GCALAB_GRAPHICS_MODE:
			break;
	}
	return;
}

/* PrintUsage(): Prints the help menu
 */
void GCALab_PrintUsage()
{
	printf("CGALab - Options");
	printf("\t [-i,--interactive]\n\t\t : start in interactive mode\n");
	printf("\t [-b,--batch]\n\t\t : start in batch mode\n");
	printf("\t [-l,--load] filename\n\t\t : load *.gca file\n");			
	printf("\t [-o,--save] filename\n\t\t : save results in *.gca file\n");
}

/* GCALab_CommandPrompt(): Prompts the user to enter a command in text mode, the
 *                         users input string is then tokenised with whitespace as the delimiter.
 */
char ** GCALab_CommandPrompt(int *numargs)
{
	char buffer[256];
	int i,j,k;
	int len;
	int argc;
	char** argv;
	char c;
	
	i = 0;
	/*show the prompt to the user*/	
	fprintf(stdout," -->");
	fflush(stdout);
	do {
		c = fgetc(stdin);
		buffer[i] = c;
		i++;
	} while(c != '\n' && i <= 255);
	len = i;
	buffer[len-1] = '\0';
	/*count spaces*/
	argc = 0;
	for (i=1;i<len;i++)
	{
		if (((buffer[i] == ' ') || (buffer[i] == '\0')) && (buffer[i-1] != ' '))		
		{
			argc++;
		}
	}
	/*allocate memory for the commands*/
	if (!(argv = (char **)malloc(argc*sizeof(char*))))
	{
		return NULL;
	}

	for (i=0;i<argc;i++)
	{
		if (!(argv[i] = (char*)malloc(128*sizeof(char))))
		{
			return NULL;
		}
	}

	j = 0;
	k = 0;
	/*parse the commandline*/
	for (i=0;i<len-1;i++)
	{
		if (buffer[i] != ' ')
		{
			argv[j][k] = buffer[i];
			k++;
			if (buffer[i+1] == ' ' || buffer[i+1] == '\0')
			{
				argv[j][k] = '\0';
				k = 0;
				j++;
			}
		}
	}
	*numargs = argc;
	return argv;
}

/* GCALab_ReadScriptCommand(): Prompts the user to enter a command in text mode, the
 *                         users input string is then tokenised with whitespace as the delimiter.
 */
char ** GCALab_ReadScriptCommand(char * filename, int *numargs)
{
	static char **cmds = NULL;
	static FILE *fp = NULL;
	static int numcmds = 0;
	static int curcmd = 0;
	char buffer[256];
	int i,j,k;
	int len;
	int argc;
	char** argv;
	char c;
	i = 0;
	if (filename != NULL)
	{
		if (!(fp = fopen(filename,"r")))
		{
			*numargs = 0;
			return NULL;
		}

		while (!feof(fp))
		{
			c = fgetc(fp);
			if (c == '\n'){
				numcmds++;
			}
		}
		rewind(fp);

		if (!(cmds = (char **)malloc(numcmds*sizeof(char*))))
		{
			*numargs = 0;
			return NULL;
		}

		for (i=0;i<numcmds;i++)
		{
			if (!(cmds[i] = (char*)malloc(256*sizeof(char))))
			{
				*numargs = 0;
				return NULL;
			}
		}

		for (i=0;i<numcmds;i++)
		{
			c = fgetc(fp);
			j = 0;
			while (c != '\n')
			{
				cmds[i][j] = c;
				c = fgetc(fp);
				j++;
			}
			cmds[i][j] = '\0';
		}
		fclose(fp);

	}
	i = 0;
	if (curcmd < numcmds)
	{

		do {
			c = cmds[curcmd][i];
			buffer[i] = c;
			i++;
		} while(c != '\0' && i <= 255);
		len = i;
		buffer[len-1] = '\0';
		/*count spaces*/
		argc = 0;
		for (i=1;i<len;i++)
		{
			if (((buffer[i] == ' ') || (buffer[i] == '\0')) && (buffer[i-1] != ' '))		
			{
				argc++;
			}
		}
		/*allocate memory for the commands*/
		if (!(argv = (char **)malloc(argc*sizeof(char*))))
		{
			return NULL;
		}

		for (i=0;i<argc;i++)
		{
			if (!(argv[i] = (char*)malloc(128*sizeof(char))))
			{	
				return NULL;
			}
		}
	
		j = 0;
		k = 0;
		/*parse the commandline*/
		for (i=0;i<len-1;i++)
		{
			if (buffer[i] != ' ')
			{	
				argv[j][k] = buffer[i];
				k++;
				if (buffer[i+1] == ' ' || buffer[i+1] == '\0')
				{
					argv[j][k] = '\0';
					k = 0;
					j++;
				}
			}
		}
		curcmd++;
		*numargs = argc;
		return argv;
	}
	else
	{
		*numargs = -1;
		return NULL;
	}
}

/* GCALab_InitCL_Options(): Initialises the options struct with defaults values
 */
void GCALab_InitCL_Options(GCALab_CL_Options* opts)
{
	opts->mode = GCALAB_DEFAULT_MODE;
	opts->load_CA = 0;
	opts->save_CA = 0;
	opts->CAInputFilename = "";
	opts->CAOutputFilename = "";
}

/* GCALab_ParseCommandLineArgs(): parses the command line arguments and returns a struct
 *                         of user options. 
 */
GCALab_CL_Options* GCALab_ParseCommandLineArgs(int argc, char **argv)
{
	int i;
	GCALab_CL_Options* CL_opt;
	
	if(!(CL_opt = (GCALab_CL_Options*)malloc(sizeof(GCALab_CL_Options))))
	{
		return NULL;
	}
	
	GCALab_InitCL_Options(CL_opt);
	
	for (i=1;i<argc;i++)
	{
		/*short format options -abcd*/
		if (argv[i][0] == '-' && argv[i][1] != '-')
		{
			int j;
			char *optslist;
			j=1;
			optslist = argv[i];
			while (optslist[j] != '\0')
			{
				
				switch (optslist[j])
				{
					case 'i':
						CL_opt->mode = GCALAB_TEXT_MODE;
						break;
					case 'b':
						CL_opt->mode = GCALAB_BATCH_MODE;
						CL_opt->ScriptFile = argv[++i];
						break;
					case 'g':
						CL_opt->mode = GCALAB_GRAPHICS_MODE;
						break;
					case 'l':
						if (optslist[j+1] == '\0')
						{
							CL_opt->load_CA = 1;
							CL_opt->CAInputFilename = argv[++i];	 
						}
						else
						{
							GCALab_PrintUsage();
							return NULL;
						}
						break;
					case 'o':
						if (optslist[j+1] == '\0')
						{
							CL_opt->save_CA = 1;
							CL_opt->CAOutputFilename = argv[++i];	
						}
						else
						{
							GCALab_PrintUsage();
							return NULL;
						}
						break;
				}
				
				j++;
			}
		} /*long format options --some_option*/
		else if(argv[i][0] == '-' && argv[i][1] == '-')
		{
			if (!strcmp(argv[i],"--interactive"))
			{
				CL_opt->mode = GCALAB_TEXT_MODE;
			}
			else if (!strcmp(argv[i],"--batch"))
			{
				CL_opt->mode = GCALAB_BATCH_MODE;
				CL_opt->ScriptFile = argv[++i];
			}
			else if(!strcmp(argv[i],"--graphics"))
			{
				CL_opt->mode = GCALAB_GRAPHICS_MODE;
			}
			else if (!strcmp(argv[i],"--load"))
			{
				CL_opt->load_CA = 1;
				CL_opt->CAInputFilename = argv[++i];	
			}
			else if(!strcmp(argv[i],"--save"))
			{
				CL_opt->save_CA = 1;
				CL_opt->CAOutputFilename = argv[++i];	
			}
			else /*unknown option*/
			{
				GCALab_PrintUsage();
				return NULL;
			}
			
		}
	}
	
	return CL_opt;	
}

/* GCALab_NewWorkSpace(): Creates a new processing workspace, the user can set the limit 
 *                        on the number of CA objects the workspace can hold.
 */
char GCALab_NewWorkSpace(int GCALimit)
{
	GCALab_WS *new_ws;
	int rc;

	if (GCALab_numWS < GCALAB_MAX_WORKSPACES)
	{
		if (!(new_ws = (GCALab_WS*)malloc(sizeof(GCALab_WS))))
		{
			return GCALAB_MEM_ERROR;
		}
		/*initialise GCA memory*/
		if(!(new_ws->GCAList = (GraphCellularAutomaton **)malloc(GCALimit*sizeof(GraphCellularAutomaton*))))
		{
			return GCALAB_MEM_ERROR;
		}
		if(!(new_ws->GCAGeometry = (mesh **)malloc(GCALimit*sizeof(mesh*))))
		{
			return GCALAB_MEM_ERROR;
		}

		new_ws->numGCA = 0;
		if (GCALimit > 0)
		{
			new_ws->maxGCA = GCALimit;
		}
		else
		{
			new_ws->maxGCA = GCALAB_DEFAULT_MAX_GCA;
		}

		/*result data memory*/
		new_ws->numresults = 0;
		if (!(new_ws->results = (GCALabOutput**)malloc(GCALAB_MAX_RESULTS*sizeof(GCALabOutput*))))
		{
			return GCALAB_MEM_ERROR;
		}
		/*set up the command queue*/
		if (!(new_ws->commandqueue = (unsigned int *)malloc(GCALAB_COMMAND_BUFFER_SIZE*sizeof(unsigned int))))
		{
			return GCALAB_MEM_ERROR;
		}
		memset((void*)(new_ws->commandqueue),0,GCALAB_COMMAND_BUFFER_SIZE*sizeof(unsigned int));
		if (!(new_ws->commandtarget = (unsigned int *)malloc(GCALAB_COMMAND_BUFFER_SIZE*sizeof(unsigned int))))
		{
			return GCALAB_MEM_ERROR;
		}
		memset((void*)(new_ws->commandtarget),0,GCALAB_COMMAND_BUFFER_SIZE*sizeof(unsigned int));
		if (!(new_ws->commandparams = (char***)malloc(GCALAB_COMMAND_BUFFER_SIZE*sizeof(char**))))
		{
			return GCALAB_MEM_ERROR;
		}
		memset((void*)(new_ws->commandparams),0,GCALAB_COMMAND_BUFFER_SIZE*sizeof(char**));
		if (!(new_ws->numparams = (int *)malloc(GCALAB_COMMAND_BUFFER_SIZE*sizeof(int))))
		{
			return GCALAB_MEM_ERROR;
		}
		memset((void*)(new_ws->numparams),0,GCALAB_COMMAND_BUFFER_SIZE*sizeof(int));


		new_ws->numcommands = 0;
		new_ws->qhead = 0;
		new_ws->qtail = 0;

		GCALab_Global[GCALab_numWS] = new_ws;
		/*init worker thread*/
		new_ws->state = GCALAB_WS_STATE_IDLE;
		/*because the main thread updates the queue*/
		pthread_mutex_init(&(new_ws->wslock),NULL);

		rc = pthread_create(&(new_ws->worker),NULL,GCALab_Worker,(void*)(unsigned long long)GCALab_numWS);
		if (rc)
		{
			return GCALAB_THREAD_ERROR;
		}

		GCALab_numWS++;
		return GCALAB_SUCCESS;
	}
	else
	{
		return GCALAB_INVALID_WS_ERROR; 
	}
}

/*  GCALab_LockWS(): aquire a lock on the given workspace
 */
void GCALab_LockWS(unsigned char ws_id)
{
	pthread_mutex_lock(&(GCALab_Global[ws_id]->wslock));
}

/*  GCALab_UnLockWS(): release a lock on the given workspace
 */
void GCALab_UnLockWS(unsigned char ws_id)
{
	pthread_mutex_unlock(&(GCALab_Global[ws_id]->wslock));
}

/* GCALab_QueueCommand(): appends the given command string to the command queue
 *                        of the given workspace.
 */
char GCALab_QueueCommand(unsigned char ws_id,unsigned int command_id,unsigned int target_id,char **params,int numparams)
{
	unsigned int ind;
	/*aquire lock on the workspace*/
	GCALab_LockWS(ws_id);
	if (WS(ws_id)->numcommands == GCALAB_COMMAND_BUFFER_SIZE)
	{
		GCALab_UnLockWS(ws_id);
		return GCALAB_QUEUE_FULL;
	}
	ind = WS(ws_id)->qtail;
	WS(ws_id)->commandqueue[ind] = command_id;
	WS(ws_id)->commandtarget[ind] = target_id;
	WS(ws_id)->commandparams[ind] = params;
	WS(ws_id)->numparams[ind] = numparams;
	WS(ws_id)->qtail = (WS(ws_id)->qtail+1)%GCALAB_COMMAND_BUFFER_SIZE;
	WS(ws_id)->numcommands++;
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}
/* GCALab_CancelCommand(): cancels the the command located at address index of 
 *                         the command queue of the given workspace.
 */
char GCALab_CancelCommand(unsigned char ws_id,unsigned int index)
{
	GCALab_LockWS(ws_id);
	if ((index >= 0) && (index < GCALAB_COMMAND_BUFFER_SIZE))
	{
		WS(ws_id)->commandqueue[index] = GCALAB_NOP;
		WS(ws_id)->commandtarget[index] = 0;
		/*clean up parameter memory if there were any*/
		if (WS(ws_id)->numparams[index] > 0)
		{
			free(WS(ws_id)->commandparams[index]);
			WS(ws_id)->commandparams[index] = NULL;
			WS(ws_id)->numparams[index] = 0;
		}

	}
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/* GCALab_ProcessCommandQueue(): Tells workspace to continue processing
 */
char GCALab_ProcessCommandQueue(unsigned char ws_id)
{
	GCALab_LockWS(ws_id);
	WS(ws_id)->state = GCALAB_WS_STATE_IDLE;
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/* GCALab_PrauseCommandQueue(): Tells workspace to halt processing
 */
char GCALab_PauseCommandQueue(unsigned char ws_id)
{
	GCALab_LockWS(ws_id);
	WS(ws_id)->state = GCALAB_WS_STATE_PAUSED;
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/* GCALab_PrintCommandQueue(): Tells workspace print all queued commands
 */
char GCALab_PrintCommandQueue(unsigned char ws_id)
{
	int i,j,k;
	GCALab_LockWS(ws_id);
	j = WS(ws_id)->qhead;
	fprintf(stdout,"%d %d\n",WS(ws_id)->qhead,WS(ws_id)->qtail);
	fprintf(stdout,"Priority\tCode\tTarget\t#Parameters\n");
	fprintf(stdout,"-----------------------------------\n");
	for(i=0;i < WS(ws_id)->numcommands;i++)
	{
		k = (i+j)%GCALAB_COMMAND_BUFFER_SIZE;
		fprintf(stdout,"\t%d\t%u\t%u\t%d\n",i,WS(ws_id)->commandqueue[k],WS(ws_id)->commandtarget[k],WS(ws_id)->numparams[k]);
	}
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/* GCALab_GetState(): gets the state flag of a workspace
 */
unsigned int GCALab_GetState(unsigned char ws_id)
{
	unsigned int state;
	GCALab_LockWS(ws_id);
	state = GCALab_Global[ws_id]->state;
	GCALab_UnLockWS(ws_id);
 	return state;
}

/* GCALab_SetState(): sets the state flag of a workspace safely
 */
void GCALab_SetState(unsigned char ws_id,unsigned int state)
{
	GCALab_LockWS(ws_id);
	GCALab_Global[ws_id]->state = state;
	GCALab_UnLockWS(ws_id);
	return;
}

/* GCALab_ShutDown(): Nicely and humainly kills the session
 */
void GCALab_ShutDown(char rc)
{
	int i;
	/*TODO: check if any workspaces are in processing state*/
	/*if so then prompt the user*/
	if (GCALab_mode == GCALAB_BATCH_MODE)
	{
		int sum;
		sum = -1;
		while (sum != 0)
		{
			sum = 0;
			for (i=0;i<GCALab_numWS;i++)
			{
				GCALab_LockWS(i);
				sum += WS(i)->numcommands;
				GCALab_UnLockWS(i);
			}
		}

	}
	for (i=0;i<GCALab_numWS;i++)
	{
		GCALab_SetState(i,GCALAB_WS_STATE_EXITING);
	}
	exit(0);
}

/* GCALab_PrintHelp(): Prints the help menu
 */
void GCALab_PrintHelp(void)
{
	int i;

	fprintf(stdout,"\nCommand List:\n");
	fprintf(stdout,"-------------\n");
	fprintf(stdout,"Command:\t[args]\tDescription\n");
	for (i=0;i<GCALab_numCmds;i++)
	{
		fprintf(stdout,"%s:\t[%s]\t%s\n",GCALab_Cmds[i].id,GCALab_Cmds[i].args,GCALab_Cmds[i].desc);
	}

	return;
}

/* GCALab_PrintOperations(): Print list of queuable compute functions
 */
void GCALab_PrintOperations(void)
{
	int i;

	fprintf(stdout,"\nOperations List:\n");
	fprintf(stdout,"-----------------\n");
	for (i=0;i<GCALab_numOps;i++)
	{
		fprintf(stdout,"Synopsis: %s %s\n",GCALab_Ops[i].id,GCALab_Ops[i].args);
		fprintf(stdout,"Description: %s\n",GCALab_Ops[i].desc);
	}
	return;
}

/* GCALab_PrintWorkSpace(): prints summart information about the 
 *                          given workspace.
 */
char GCALab_PrintWorkSpace(unsigned char ws_id)
{
	int i;

	fprintf(stdout,"Work Space ID: %hhu\n",ws_id);
	fprintf(stdout,"GCA (%d): ",GCALab_Global[ws_id]->numGCA);
	for (i=0;i<GCALab_Global[ws_id]->numGCA;i++)
	{
		fprintf(stdout,"(%d,%u,%u) ",i,GCALab_Global[ws_id]->GCAList[i]->params->rule,GCALab_Global[ws_id]->GCAList[i]->params->N);
	}
	fprintf(stdout,"\n");
	
	fprintf(stdout,"Results (%d): ",GCALab_Global[ws_id]->numresults);
	for (i=0;i<GCALab_Global[ws_id]->numresults;i++)
	{
		fprintf(stdout,"(%d,%d,%s)",i,GCALab_Global[ws_id]->results[i]->type,GCALab_Global[ws_id]->results[i]->id);
	}
	fprintf(stdout,"\n");

	fprintf(stdout,"Queued Commands (%d):\n",GCALab_Global[ws_id]->numcommands);
	fprintf(stdout,"Current State: ");
	fprintf(stdout,"%s\n",statenames[GCALab_Global[ws_id]->state]);
	return GCALAB_SUCCESS;
}

/* GCALab_ListWorkSpaces(): Lists summary info about all work spaces
 */
char GCALab_ListWorkSpaces(void)
{
	int i;
	fprintf(stdout,"ID\tCA\tQL\tST\n");
	for (i=0;i<GCALab_numWS;i++)
	{
		fprintf(stdout,"%d\t%d\t%u\t",i,GCALab_Global[i]->numGCA,GCALab_Global[i]->numcommands);
		fprintf(stdout,"%c\n",stateInitials[GCALab_Global[i]->state]);
	}
	return GCALAB_SUCCESS;
}

/* GCALab_PrintCA(): Prints Details of Selected CA including the current configuration
 */
char GCALab_PrintCA(unsigned char ws_id,unsigned int gca_id)
{
	int i;
	GraphCellularAutomaton *GCA;
	if (gca_id < 0 || gca_id >= WS(ws_id)->numGCA)
	{
		return GCALAB_INVALID_OPTION;
	}

	GCA = WS(ws_id)->GCAList[gca_id];

	fprintf(stdout,"#cells: %u\n",GCA->params->N);
	fprintf(stdout,"k-neighbourhood: %u\n",GCA->params->k);
	fprintf(stdout,"rule code: %u\n",GCA->params->rule);
	fprintf(stdout,"rule type: %u\n",GCA->params->rule_type);
	fprintf(stdout,"time-window: %d\n",GCA->params->WSIZE);
	fprintf(stdout,"Configuration at t = 0\n");
	for(i=0;i<GCA->params->N;i++)
	{
		fprintf(stdout,"%c",cellsymbols[GetCellStatePacked_external(GCA,GCA->ic,i)]);
	}
	fprintf(stdout,"\n");
	fprintf(stdout,"Configuration at t = %d\n",GCA->t);
	for(i=0;i<GCA->params->N;i++)
	{
		fprintf(stdout,"%c",cellsymbols[GetCellStatePacked_external(GCA,GCA->config,i)]);
	}
	fprintf(stdout,"\n");

	return GCALAB_SUCCESS;
}

/* GCALab_PrintSTP(): Prints the space-time pattern of the select CA
 */
char GCALab_PrintSTP(unsigned char ws_id,unsigned int gca_id)
{
	int i,l,j;
	GraphCellularAutomaton *GCA;

	if (gca_id < 0 || gca_id >= WS(ws_id)->numGCA)
	{
		return GCALAB_INVALID_OPTION;
	}
	

	GCA = WS(ws_id)->GCAList[gca_id];
	l = (GCA->t >= GCA->params->WSIZE) ? GCA->params->WSIZE - 1 :  GCA->t;
	for (i=l;i>0;i--)
	{
		for (j=0;j<GCA->params->N;j++)
		{
			fprintf(stdout,"%c",cellsymbols[GetCellStatePacked(GCA,j,i)]);	
		}
		fprintf(stdout,"\n");
	}

	return GCALAB_SUCCESS;
}

/* GCALab_PrintResults(): Prints compute results
 */
char GCALab_PrintResults(unsigned char ws_id,unsigned int res_id)
{
	GCALabOutput *res;
	int i;
	if (res_id < 0 || res_id >= WS(ws_id)->numresults)
	{
		return	GCALAB_INVALID_OPTION;
	}
	
	res = WS(ws_id)->results[res_id];
	fprintf(stdout,"type: %d\n",res->type);
	fprintf(stdout,"id: %s\n",res->id);
	fprintf(stdout,"data length: %u\n",res->datalen);
	switch(res->type)
	{
		case FLOAT32:
		{	
			float *d;
			d = (float*)res->data;
			for (i=0;i<res->datalen;i++)
			{
				fprintf(stdout,"%f\n",d[i]);
			}
		}
			break;
		case UINT32:
		{	
			unsigned int *d;
			d = (unsigned int*)res->data;
			for (i=0;i<res->datalen;i++)
			{
				fprintf(stdout,"%u\n",d[i]);
			}
		}
			break;
		case CHUNK:
		{	
			chunk *d;
			d = (chunk*)res->data;
			for (i=0;i<res->datalen;i++)
			{
				fprintf(stdout,"%X\n",d[i]);
			}
		}
			break;
	}
	return GCALAB_SUCCESS;
}

/* strvncpy(): makes a copy of an array of strings 
 */
char** strvncpy(char **strv,int c, int n)
{
	char** strv_cp;
	int i;

	if(!(strv_cp = (char**)malloc(c*sizeof(char*))))
	{
		return NULL;
	}

	for (i=0;i<c;i++)
	{
		if(!(strv_cp[i] = (char*)malloc(n*sizeof(char))))
		{
			return NULL;
		}
	}

	for (i=0;i<c;i++)
	{
		strncpy(strv_cp[i],strv[i],n);
	}

	return strv_cp;
}

#ifdef WITH_GRAPHICS
/* GCALab_Graphics_Init(): Initialises the graphics settings
 */
void GCALab_Graphics_Init(void)
{
	/*background colour*/
	glClearColor(0.5,0.5,0.5,1.0);
	/*enable depth test, blending, back face culling, lighting, and anti-aliasing*/
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_MULTISAMPLE); 
	/*lighting model*/
	glEnable(GL_LIGHT0); 	
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_POSITION, position)   ; 
	/*smooth shading model*/
	glShadeModel(GL_SMOOTH);
	/*fog effect stuff*/
	glFogi(GL_FOG_MODE,GL_LINEAR);
	glFogfv(GL_FOG_COLOR,fogCol);
	glFogf(GL_FOG_DENSITY,0.55);
	glHint(GL_FOG_HINT,GL_NICEST);
	glFogf(GL_FOG_START,1.0);
	glFogf(GL_FOG_END,50.0);
	glEnable(GL_FOG);
	quad = gluNewQuadric();
	return;
}

/* Render the GCA scene
 */
void GCALab_Graphics_Display(void)
{
	/*clear screen for a new render*/
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	
	/*transformations*/
	glTranslatef(translateX,translateY,zoom-10.5f);
	glRotatef(phi+30.0f,1.0f,0.0f,0.0f);
	glRotatef(theta,0.0f,1.0f,0.0f);
	
	/* draw the reflected scene first*/
	glLightfv(GL_LIGHT0,GL_POSITION,positionMirror);
	glPushMatrix();
	glFrontFace(GL_CW);
	glScalef(1.0,-1.0,1.0);
	glPushMatrix();
	glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
	GCALab_Graphics_DrawBBox();
	glPopMatrix();
	if (GCALab_numWS > 0)
	{
		if (WS(cur_ws)->numGCA > 0)
		{
			glPushMatrix();
			glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
			GCALab_Graphics_DrawGCA(cur_ws,cur_gca);
			glPopMatrix();
		}
	} 
	glFrontFace(GL_CCW);
	glPopMatrix();
    

	GCALab_Graphics_DrawGrid();

	/* draw the "real" scene... but how do we define real?*/
	glLightfv(GL_LIGHT0, GL_POSITION, position); 
	
	glPushMatrix();
	glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
	GCALab_Graphics_DrawBBox();
	glPopMatrix();
	if (GCALab_numWS > 0)
	{
		if (WS(cur_ws)->numGCA > 0)
		{
			glPushMatrix();
			glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
			GCALab_Graphics_DrawGCA(cur_ws,cur_gca);
			glPopMatrix();
		}
	} 

	glutSwapBuffers();
}

/* GCALab_Graphics_DrawGCA(): draws the selected graph cellular automaton
 *                            to the OpenGL scene.
 */
void GCALab_Graphics_DrawGCA(unsigned int ws_id,unsigned int gca_id)
{
	float *n;
	mesh *m;
	int i;
	GraphCellularAutomaton *gca;
	m = WS(ws_id)->GCAGeometry[gca_id];
	gca = WS(ws_id)->GCAList[gca_id];
	
	/*we assume that we are dealing with a triangular mesh*/
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	/*TODO: expand this to support other meshes*/
	glBegin(GL_TRIANGLES);

	for (i=0;i<m->fList->numFaces;i++)
	{
		/*change material properties based on cell state*/
		float * col;
		if (lutColmode)
		{
			col = lutColf[GetNeighbourhood_config(gca,i,0)];
		}
		else
		{
			col = cellColf[GetCellStatePacked(gca,i,0)];
		}
		GLfloat ambdif[4] = {col[0],col[1],col[2],1.0};
		GLfloat em[4] = {0.0,0.0,0.0,1.0};
		GLfloat spec[4] = {0.7,0.7,0.7,0.5};
		GLfloat shininess  = 128.0;
		glMaterialfv(GL_FRONT,GL_EMISSION,em);   
		glMaterialfv(GL_FRONT,GL_SPECULAR, spec);
		glMaterialf(GL_FRONT,GL_SHININESS,shininess);
		glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,ambdif);    
		/*compute the cell normal*/
		n = Normal_f(m->vList->verts + 3*(m->fList->faces[3*i]),m->vList->verts + 3*(m->fList->faces[3*i+1]),m->vList->verts + 3*(m->fList->faces[3*i+2]));
		/*draw the cell*/
		glNormal3fv(n);
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i]));
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+1]));
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+2]));
		free(n);
	}
	glEnd();
}

/* GCALab_Graphics_DrawBBox(): draws a bounding box
 */
void GCALab_Graphics_DrawBBox(void)
{
    glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor4f(1,0,0,1);
	glVertex3f(0,0,0);
	glVertex3f(1,0,0);
	
	glColor4f(0,1,0,1);
	glVertex3f(0,0,0);
	glVertex3f(0,1,0);

	glColor4f(0,0,1,1);
	glVertex3f(0,0,0);
	glVertex3f(0,0,1);
	glEnd();
    glEnable(GL_LIGHTING);
}

/* GCALab_Graphics_DrawGrid(): draws a grid to orientate the user
 */
void GCALab_Graphics_DrawGrid(void)
{
	int lim,i,j;
	/* draw the transparent surface*/
    glDisable(GL_LIGHTING);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	lim = 200;

	glColor4f(0.0,0.0,0.0,0.7);
	glBegin(GL_QUADS);
		glVertex3f((float)(-lim/2),0.0,(float)(-lim/2));
		glVertex3f((float)(-lim/2),0.0,(float)(lim/2));
		glVertex3f((float)(lim/2),0.0,(float)(lim/2));
		glVertex3f((float)(lim/2),0.0,(float)(-lim/2));
	glEnd();

	/* draw some gridlines*/
	glColor4f(0.0,0.0,0.0,0.9);
	for (i=0;i<lim;i++)
	{
		for (j=0;j<lim;j++)
		{
			glBegin(GL_LINE_LOOP);
				glVertex3f((float)(-lim/2+i),0.2,(float)(-lim/2+j));
				glVertex3f((float)(-lim/2+i),0.2,(float)(-lim/2+1+j));
				glVertex3f((float)(-lim/2+1+i),0.2,(float)(-lim/2+1+j));
				glVertex3f((float)(-lim/2+1+i),0.2,(float)(-lim/2+j));
			glEnd();		
		}
	}
	glEnable(GL_LIGHTING);
}

/* GCALab_Graphics_keyPressed():handles standard key press events
 */
void GCALab_Graphics_KeyPressed (unsigned char key, int x, int y)
{
	char rc;
	switch(key)
	{
		case 'c': /*enter a command directly*/
		{
			char **userinput;
			int numargs;
			char* cmd;
			int i;
			/*get user input*/
			userinput = GCALab_CommandPrompt(&numargs);
			/*ensure input is valid*/
			rc = GCALab_TestPointer((void*)userinput);
			GCALab_HandleErr(rc);
			
			cmd = userinput[0];
			for (i=0;i<GCALab_numCmds;i++)
			{
				if (!strcmp(cmd,GCALab_Cmds[i].id))
				{
					rc = (*(GCALab_Cmds[i].f))(numargs,userinput);
					GCALab_HandleErr(rc);
				}
			}
			/*clean up user args*/
			for (i=0;i<numargs;i++)
			{
				free(userinput[i]);
			}
			free(userinput);
		}
			break;
		case '>':
		{
			unsigned int new_ws;
			/*next workspace*/
			new_ws = cur_ws+1;
			rc = GCALab_ValidWSId(new_ws);
			/*only update if the id was valid*/
			if (rc == GCALAB_INVALID_WS_ERROR) return;
			cur_ws = new_ws;
			fprintf(stdout,"Current Workspace is ID = %d\n",cur_ws);
		}
			break;
		case '<':
		{
			unsigned int new_ws;
			/*prev workspace*/
			new_ws = cur_ws-1;
			rc = GCALab_ValidWSId(new_ws);
			/*only update if the id was valid*/
			if (rc == GCALAB_INVALID_WS_ERROR) return;
			cur_ws = new_ws;
			fprintf(stdout,"Current Workspace is ID = %d\n",cur_ws);
		}	
			break;
		case ']':
			/*next gca*/
			cur_gca = (cur_gca >= WS(cur_ws)->numGCA-1) ?  WS(cur_ws)->numGCA-1  : cur_gca + 1;
			fprintf(stdout,"GCA %d selected\n",cur_gca);
			break;
		case '[':
			/*prev gca*/
			cur_gca = (cur_gca <= 0) ?  0 : cur_gca - 1;
			fprintf(stdout,"GCA %d selected\n",cur_gca);
			break;
		case '}':
			/*next result*/
			cur_res = (cur_res >= WS(cur_ws)->numresults-1) ?  WS(cur_ws)->numresults-1 : cur_res + 1;
			fprintf(stdout,"Result %d selected\n",cur_res);
			break;
		case '{':
			/*prev result*/
			cur_res = (cur_res <= 0) ? 0 : cur_res - 1;
			fprintf(stdout,"Result %d selected\n",cur_res);
			break;
		case 'Q':
			/*quit GCALab*/
			GCALab_ShutDown(0); 
			break;
		case 's':
			CANextStep(WS(cur_ws)->GCAList[cur_gca]);
			break;
		case 'l':
			lutColmode = !lutColmode;
			break;
	}
}

/* GCALab_Graphics_SpecialKeyPressed(): handles special key press events
 */
void GCALab_Graphics_SpecialKeyPressed(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_LEFT: 
			theta -= 0.5;
			theta = (theta < 0.0) ? 360.0 : theta;
			break;
		case GLUT_KEY_UP: 
			zoom -= 0.1;
			break;
		case GLUT_KEY_RIGHT: 
			theta += 0.5;
			theta = (theta > 360.0) ? 0.0 : theta;
			break;
		case GLUT_KEY_DOWN: 
			zoom += 0.1;
			break;
		case GLUT_KEY_HOME:
			translateY += 0.1;
			break;
		case GLUT_KEY_END:
			translateY -= 0.1;
			break;
		case GLUT_KEY_F12:
			translateX += 0.1;
			break;
		case GLUT_KEY_INSERT:
			translateX -= 0.1;
			break;
		case GLUT_KEY_PAGE_UP:
			phi+=0.5;
			break;
		case GLUT_KEY_PAGE_DOWN:
			phi-=0.5;
			break;
	}
}
/* GCALab_Graphics_Mouse(): Handles mouse events
 */
void GCALab_Graphics_MouseClick(int button,int state,int x, int y)
{
	switch(button)
	{
		case GLUT_LEFT_BUTTON:
			switch(state)
			{
				case GLUT_DOWN:
					moving = GCALAB_GRAPHICS_INTERACTION_ROTATE;
					xprev = x;
					yprev = y;
					break;
				case GLUT_UP:
					moving = GCALAB_GRAPHICS_INTERACTION_NONE;
					break;
			}
			break;
		case GLUT_MIDDLE_BUTTON:
			switch(state)
			{
				case GLUT_DOWN:
					moving = GCALAB_GRAPHICS_INTERACTION_ZOOM;
					xprev = x;
					yprev = y;
					break;
				case GLUT_UP:
					moving = GCALAB_GRAPHICS_INTERACTION_NONE;
					break;
			}
			break;
		case GLUT_RIGHT_BUTTON:
			switch(state)
			{
				case GLUT_DOWN:
					moving = GCALAB_GRAPHICS_INTERACTION_TRANSLATE;
					xprev = x;
					yprev = y;
					break;
				case GLUT_UP:
					moving = GCALAB_GRAPHICS_INTERACTION_NONE;
					break;
			}
			break;		
	}
}

/* GCALab_Graphics_MouseMove(): handles mouse movement
 */
void GCALab_Graphics_MouseMove(int x, int y)
{
	switch (moving)
	{
		case GCALAB_GRAPHICS_INTERACTION_NONE:
			break;
		case GCALAB_GRAPHICS_INTERACTION_ROTATE:
			theta += ((float)(x - xprev))*0.1;
			phi += ((float)(y - yprev))*0.1;
			phi = (phi < -25.0) ? -25.0 : ((phi > 60.0) ? 60.0 : phi);
			xprev = x;
			yprev = y;
			break;
		case GCALAB_GRAPHICS_INTERACTION_ZOOM:
			zoom += ((float)(y - yprev))*0.05;
			xprev = x;
			yprev = y;
			break;
		case GCALAB_GRAPHICS_INTERACTION_TRANSLATE:
			translateX -= ((float)(x - xprev))*0.05;
			translateY += ((float)(y - yprev))*0.05;
			xprev = x;
			yprev = y;
			break;
	}
}

/* GCALab_Graphics_Reshape(): Sets up projection matrix for 
 *                            perspective projection
 */
void GCALab_Graphics_Reshape(int w, int h)
{
	if (h == 0)
	{
		h = 1;
	}
	glViewport(0,0,w,h);
	/* define projection matirx*/
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f,((float)w)/((float)h),0.1f,1000000.0f);
	/* define transformation matrix from model to world coordinates*/
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(0.0,0.0,-10.0,0.0,0.0,0.0,0.0,1.0,0.0);
	return;
}
/* GCALab_Graphics_Timer(): timer callback called roughly 60
 *                          times per second.
 */
void GCALab_Graphics_Timer(int x)
{
	glutPostRedisplay();
	glutTimerFunc(1000.0/60.0,GCALab_Graphics_Timer,x);
}
#endif
/*menu commands*/

/* GCALab_CMD_NewWorkSpace(): GCALab command to create a new
 *                            workspace
 */
char GCALab_CMD_NewWorkSpace(int argc, char **argv)
{
	int lim;
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		lim = atoi(argv[1]);
		rc = GCALab_NewWorkSpace(lim);
		if (rc != GCALAB_SUCCESS) return rc;
		cur_ws = GCALab_numWS - 1;
		printf("New Workspace created! ID = %d\n",cur_ws);
		return GCALAB_SUCCESS;
	}
}

/* GCALab_CMD_PrintWorkSpace(): GCALab command to print
 *                              the current workspace
 */
char GCALab_CMD_PrintWorkSpace(int argc, char **argv)
{
	return GCALab_PrintWorkSpace(cur_ws);
}

/* GCALab_CMD_ListWorkSpace(): GCALab command to print summary 
 *                            information on all workspaces
 */
char GCALab_CMD_ListWorkSpaces(int argc, char **argv)
{
	return GCALab_ListWorkSpaces();
}

/* GCALab_CMD_PrintHelp(): GCALab command to print help menu
 */
char GCALab_CMD_PrintHelp(int argc, char **argv)
{
	GCALab_PrintHelp();
	return GCALAB_SUCCESS;
}

/* GCALab_CMD_PrintOperations(): GCALab command to print function
 *                               list
 */
char GCALab_CMD_PrintOperations(int argc, char **argv)
{
	GCALab_PrintOperations();
	return GCALAB_SUCCESS;
}

/* GCALab_CMD_ChangeWorkSpace(): GCALab command to set focus to
 *                               a given workspace
 */
char GCALab_CMD_ChangeWorkSpace(int argc, char **argv)
{
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int new_ws;
		new_ws = (unsigned int)atoi(argv[1]);
		rc = GCALab_ValidWSId(new_ws);
		/*only update if the id was valid*/
		if (rc == GCALAB_INVALID_WS_ERROR) return rc;
		cur_ws = new_ws;
		fprintf(stdout,"Current Workspace is ID = %d\n",cur_ws);
		return GCALAB_SUCCESS;
	}
}

/* GCALab_CMD_QueueCommand(): GCALab command to enqueue a 
 *                            compute operation
 */
char GCALab_CMD_QueueCommand(int argc, char **argv)
{
	char rc;
	if (argc < 3)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int cmd_code;
		unsigned int target;
		int numparams;
		char ** params;
		int i;
		/*first get the id map for the command*/
		cmd_code = GCALab_GetCommandCode(argv[1]);
		/*the target id is next*/
		target = (unsigned int)atoi(argv[2]);
			
		/*everything else is specific to the command*/
		numparams = argc - 3;
		/*make a copy*/
		params = strvncpy(argv+3,numparams,GCALAB_MAX_STRLEN);
		rc = GCALab_TestPointer((void*)params);
		if (rc == GCALAB_MEM_ERROR) return rc;
			
		/*push it to the queue*/
		return GCALab_QueueCommand(cur_ws,cmd_code,target,params,numparams);	
	}
}

/* GCALab_CMD_DeleteCommand(): GCALab command to delete a compute
 *                             operation.
 */
char GCALab_CMD_DeleteCommand(int argc, char **argv)
{
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int ind;
		ind = (unsigned int)atoi(argv[1]);
		return GCALab_CancelCommand(cur_ws,ind);
	}
}

/* GCALab_CMD_ExecuteQueue(): GCALab Command to signal a thread to 
 *                            process a command queue
 */
char GCALab_CMD_ExecuteQueue(int argc, char **argv)
{
	return GCALab_ProcessCommandQueue(cur_ws);
}

/* GCALab_CMD_StopQueue(): GCALab command to signal a thread to
 *                         halt processing of a command queue.
 */
char GCALab_CMD_StopQueue(int argc, char **argv)
{
	return GCALab_PauseCommandQueue(cur_ws);
}

/* GCALab_CMD_PrintCA(): GCALab command to print CA summary
 */
char GCALab_CMD_PrintCA(int argc, char **argv)
{
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		return GCALab_PrintCA(cur_ws,id);
	}
}

/* GCALab_CMD_PrintSTP(): GCALab command to print CA evolution
 */
char GCALab_CMD_PrintSTP(int argc, char ** argv)
{
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		return GCALab_PrintSTP(cur_ws,id);
	}
}

/* GCALab_CMD_PrintResults(); GCALab command to print result data.
 */
char GCALab_CMD_PrintResults(int argc, char **argv)
{
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		return GCALab_PrintResults(cur_ws,id);
	}
}

/* GCALab_CMD_Quit(): GCALab command to quit 
 */
char GCALab_CMD_Quit(int argc, char **argv)
{
	GCALab_ShutDown(GCALAB_SUCCESS);
}


/*compute operations*/

/* GCALab_OP_NOP(): No Operation
 */
char GCALab_OP_NOP(unsigned char ws_id,unsigned int trgt,int argc, char ** argv,GCALabOutput **res)
{
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Load(): Load a CA from file.
 */
char GCALab_OP_Load(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	char *filename;
	int i;
	GraphCellularAutomaton *GCA;
	mesh *m;
	filename = NULL;
	for (i=0;i<nparams;i++)
	{
		if(!strcmp(params[i],"-f"))
		{
			filename = params[++i];
		}
	}

	GCA = NULL;
	m = NULL;
	if (filename != NULL)
	{
		GCALab_fio_loadCA(filename,&GCA,&m);
	}
	WS(ws_id)->GCAList[WS(ws_id)->numGCA] = GCA;
	WS(ws_id)->GCAGeometry[WS(ws_id)->numGCA] = m;
	WS(ws_id)->numGCA++;
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Save(): Save CA or Result data to file
 */
char GCALab_OP_Save(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	char *filename;
	unsigned char gca_res_flag;
	int i;
	GraphCellularAutomaton *GCA;
	mesh *m;
	GCALabOutput *data;
	filename = NULL;
	gca_res_flag = 0;
	for (i=0;i<nparams;i++)
	{
		if(!strcmp(params[i],"-f"))
		{
			filename = params[++i];
		}
		else if (!strcmp(params[i],"-r"))
		{
			gca_res_flag = 0;
		}
		else if (!strcmp(params[i],"-g"))
		{
			gca_res_flag = 1;
		}
	}
	if (gca_res_flag)
	{
		GCA = WS(ws_id)->GCAList[trgt_id];
		m = WS(ws_id)->GCAGeometry[trgt_id];
		GCALab_fio_saveCA(filename,GCA,m);	
	}
	else
	{
		data = WS(ws_id)->results[trgt_id];
		GCALab_fio_saveData(filename,data->id,data->data,data->datalen,data->type);
	}
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Simulate(): direct simulation of CA evolution
 */
char GCALab_OP_Simulate(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int Tfinal;
	unsigned char reInit,ic_type;
	int i;
	char *ic_filename;
	GraphCellularAutomaton *GCA;
	reInit = 0;
	for (i=0;i<nparams;i++)
	{
		if(!strcmp(params[i],"-t"))
		{
			Tfinal = (unsigned int)atoi(params[++i]);
		}
		else if (!strcmp(params[i],"-I"))
		{
			reInit = 1;
		}
		else if (!strcmp(params[i],"-f"))
		{
			ic_filename = params[++i];
		}
		else if (!strcmp(params[i],"-c"))
		{
			ic_type = (unsigned char)atoi(params[++i]);
		}
	}
	
	GCA = WS(ws_id)->GCAList[trgt_id];
	if (reInit)
	{
		ResetCA(GCA);
		SetCAIC(GCA,NULL,ic_type);
	}
	
	CASimTSteps(GCA,Tfinal);
	return GCALAB_SUCCESS;
}

/* GCALab_OP_GCA(): Create a Graph Cellular Automaton
 */
char GCALab_OP_GCA(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	char *meshfile;
	unsigned int NCell,genus,windowsize,r;
	unsigned char r_type;
	unsigned char s,k,eca,ic_type;
	int  i;
	char rc;
	GraphCellularAutomaton *GCA;
	CellularAutomatonParameters *CAparams;
	mesh *m;

	windowsize = 0;
	meshfile = NULL;
	eca = 0;
	for (i=0;i<nparams;i++)
	{
		if(!strcmp(params[i],"-m"))
		{
			meshfile = params[++i];
		}
		else if(!strcmp(params[i],"-t"))
		{
			NCell = (unsigned int)atoi(params[++i]);
			genus = (unsigned int)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-r"))
		{
			r_type = (unsigned char)atoi(params[++i]);
			r = (unsigned int)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-s"))
		{
			s = (unsigned char)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-w"))
		{
			windowsize = (unsigned int)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-eca"))
		{
			NCell = (unsigned int)atoi(params[++i]);
			k = (unsigned char)atoi(params[++i]);
			r = (unsigned int)atoi(params[++i]);
			eca = 1;
		}
		else if (!strcmp(params[i],"-c"))
		{
			ic_type = (unsigned char)atoi(params[++i]);
		}
	}

	if (eca)
	{
		/*create an elementary CA with N cells, k-neighbourhood and wolfram code r*/
		GCA = CreateECA(NCell,k,r,windowsize);
		m = NULL;
		rc = GCALab_TestPointer((void*)GCA);
		if (rc <= 0)
		{
			return rc;
		}
	}
	else
	{
		if (meshfile != NULL)
		{
			/*this is not ideal by is ok for now*/
			m = LoadMesh(meshfile,OFF_FORMAT);
		}
		else
		{
			//printf("get here!\n");
			m = CreateMeshTopology(NCell,genus);
		}
		rc = GCALab_TestPointer((void*)m);
		if (rc <= 0)
		{
			return rc;
		}
		CAparams = CreateCAParams(m,s,r_type,r,windowsize);
		rc = GCALab_TestPointer((void*)CAparams);
		if (rc <= 0)
		{
			return rc;
		}
		GCA = CreateGCA(CAparams);
		rc = GCALab_TestPointer((void*)GCA);
		if (rc <= 0)
		{
			return rc;
		}
	}
		
	ResetCA(GCA);
	SetCAIC(GCA,NULL,ic_type);
	WS(ws_id)->GCAList[WS(ws_id)->numGCA] = GCA;
	WS(ws_id)->GCAGeometry[WS(ws_id)->numGCA] = m;
	WS(ws_id)->numGCA++;
	
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Entropy(): compute entropy meansures on a given CA
 */
char GCALab_OP_Entropy(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int numSamples;
	unsigned int T;
	unsigned type;
	int i;
	char rc;
	GraphCellularAutomaton *GCA;
	
	float *p,*logs_p, *S_i,*pt,*logs_pt,*IE,*logQ;
	unsigned char *flags;
	unsigned int *count,*countt,*wl,*Q;
	
	float S_mu,W_mu,I_mu,I_sigma;
	float *result_data;
	numSamples = 1;
	type = GCALAB_SHANNON_ENTROPY;
	for (i=0;i<nparams;i++)
	{
		if (!strcmp(params[i],"-n"))
		{
			numSamples = (unsigned int)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-t"))
		{
			T = (unsigned int)atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-e"))
		{
			type = (unsigned int)atoi(params[++i]);
		}
	}

	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
	(*res) = (GCALabOutput*)malloc(sizeof(GCALabOutput)); 
	switch(type)
	{
		case GCALAB_SHANNON_ENTROPY:
			/*allocate memory first reduce malloc calls*/
			p = (float*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(float));
			logs_p = (float*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(float));
			S_i = (float*)malloc((GCA->params->N)*sizeof(float));
			flags = (unsigned char*)malloc((GCA->params->N)*sizeof(unsigned char));
			count = (unsigned int*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(unsigned int));

			rc = GCALab_TestPointer((void*)p);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logs_p);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)S_i);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)flags);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)count);
			if (rc <= 0)
			{
				return rc;
			}
			
			/*compute avg Shannon entropy*/
			S_mu = 0.0;
			for (i=0;i<numSamples;i++)
			{
				ResetCA(GCA);
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
			}

			S_mu = S_mu/((float)numSamples);
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):S",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = S_mu;
			(*res)->data = (void*)result_data;

			/*clean up*/	
			free(p);
			free(logs_p);
			free(S_i);
			free(flags);
			free(count);
			break;
		case GCALAB_WORD_ENTROPY:
			/*allocate memory first reduce malloc calls*/
			pt = (float*)malloc((GCA->params->N)*T*sizeof(float));
			logs_pt = (float*)malloc((GCA->params->N)*T*sizeof(float));
			S_i = (float*)malloc((GCA->params->N)*sizeof(float));
			flags = (unsigned char*)malloc((GCA->params->N)*sizeof(unsigned char));
			countt = (unsigned int*)malloc((GCA->params->N)*T*sizeof(unsigned int));
			wl = (unsigned int*)malloc((GCA->params->N)*sizeof(unsigned int));
			rc = GCALab_TestPointer((void*)pt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logs_pt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)S_i);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)flags);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)countt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)wl);
			if (rc <= 0)
			{
				return rc;
			}

			/*compute avg word entropy*/
			W_mu = 0.0;
			for (i=0;i<numSamples;i++)
			{
				ResetCA(GCA);
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			    W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			}

			W_mu = W_mu/((float)numSamples);
				
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):W",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = W_mu;
			(*res)->data = (void*)result_data;
			/*clean up*/	
			free(pt);
			free(logs_pt);
			free(S_i);
			free(flags);
			free(countt);
			free(wl);
			break;
		case GCALAB_INPUT_ENTROPY:
			/*allocate memory first reduce malloc calls*/
			Q =  (unsigned int *)malloc((GCA->LUT_size)*sizeof(unsigned int));
			logQ = (float *)malloc((GCA->LUT_size)*sizeof(float));
			IE = (float *) malloc(T*sizeof(unsigned int));

			rc = GCALab_TestPointer((void*)Q);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logQ);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)IE);
			if (rc <= 0)
			{
				return rc;
			}

			/*compute avg and varience of I*/
			ResetCA(GCA);
			SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			InputEntropy(GCA,T,&I_mu,&I_sigma,Q,logQ,IE);	
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):I",trgt_id);
			(*res)->datalen = T+2;
			result_data = (float*)malloc((T+2)*sizeof(float));
			result_data[0] = I_mu;
			result_data[1] = I_sigma;
			for (i=0;i<T;i++)
			{
				result_data[i+2] = IE[i];
			}
			(*res)->data = (void*)result_data;
			/*clean up*/	
			free(Q);
			free(logQ);
			free(IE);
			break;
		case GCALAB_ALL_ENTROPY:
			/*allocate memory first reduce malloc calls*/
			p = (float*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(float));
			logs_p = (float*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(float));
			S_i = (float*)malloc((GCA->params->N)*sizeof(float));
			flags = (unsigned char*)malloc((GCA->params->N)*sizeof(unsigned char));
			count = (unsigned int*)malloc((GCA->params->N)*((unsigned int)GCA->params->s)*sizeof(unsigned int));
			pt = (float*)malloc((GCA->params->N)*T*sizeof(float));
			logs_pt = (float*)malloc((GCA->params->N)*T*sizeof(float));
			countt = (unsigned int*)malloc((GCA->params->N)*T*sizeof(unsigned int));
			wl = (unsigned int*)malloc((GCA->params->N)*sizeof(unsigned int));
			Q =  (unsigned int *)malloc((GCA->LUT_size)*sizeof(unsigned int));
			logQ = (float *)malloc((GCA->LUT_size)*sizeof(float));
			IE = (float *) malloc(T*sizeof(unsigned int));
			rc = GCALab_TestPointer((void*)p);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logs_p);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)S_i);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)flags);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)count);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)pt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logs_pt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)countt);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)wl);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)Q);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)logQ);
			if (rc <= 0)
			{
				return rc;
			}
			rc = GCALab_TestPointer((void*)IE);
			if (rc <= 0)
			{
				return rc;
			}
			
			/*compute avg Shannon entropy*/
			S_mu = 0.0;
			for (i=0;i<numSamples;i++)
			{
				ResetCA(GCA);
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
			}

			S_mu = S_mu/((float)numSamples);
					
			/*compute avg word entropy*/
			W_mu = 0.0;
			for (i=0;i<numSamples;i++)
			{
				ResetCA(GCA);
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			}

			W_mu = W_mu/((float)numSamples);
					
			/*compute avg and varience of I*/
			ResetCA(GCA);
			SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			InputEntropy(GCA,T,&I_mu,&I_sigma,Q,logQ,IE);	
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):A",trgt_id);
			(*res)->datalen = T+4;
			result_data = (float*)malloc((T+4)*sizeof(float));
			result_data[0] = S_mu;
			result_data[1] = W_mu;
			result_data[2] = I_mu;
			result_data[3] = I_sigma;
			for (i=0;i<T;i++)
			{
				result_data[i+4] = IE[i];
			}
			(*res)->data = (void*)result_data;
			/*clean up*/	
			free(p);
			free(logs_p);
			free(S_i);
			free(flags);
			free(count);
			free(pt);
			free(logs_pt);
			free(countt);
			free(wl);
			free(Q);
			free(logQ);
			free(IE);
			break;
	}

	return GCALAB_SUCCESS;
}

/* GCALab_OP_Param(): compute various CA parameters (e.g. Lambda, or Z)
 */
char GCALab_OP_Param(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int type;
	chunk range[2];
	float lambdap,Zp,Gp;
	int i;
	unsigned int samples;
	GraphCellularAutomaton *GCA;
	samples = 0;
	for (i=0;i<nparams;i++)
	{
		if (!strcmp(params[i],"-p"))
		{
			type = (unsigned int )atoi(params[++i]);
		}
		else if(!strcmp(params[i],"-l"))
		{
			range[0] = (chunk)atoi(params[++i]);
			range[1] = (chunk)atoi(params[++i]);
		}
		else if (!strcmp(params[i],"-n"))
		{
			samples = (unsigned int)atoi(params[++i]);
		}
	}

	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
	(*res) = (GCALabOutput*)malloc(sizeof(GCALabOutput)); 
			
	switch(type)
	{
		case GCALAB_LAMBDA_PARAM:
		{
			float *result_data;
			lambdap = lambda_param(GCA);
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):L",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = lambdap;
			(*res)->data = (void*)result_data;
		}	
			break;
		case GCALAB_Z_PARAM:
		{
			float *result_data;
			Zp = Z_param(GCA);
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):Z",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = Zp;
			(*res)->data = (void*)result_data;
		}
			break;
		case GCALAB_G_PARAM:
		{
			float *result_data;
			if (samples > 0)
			{
				Gp = G_density(GCA,NULL,samples);
			}
			else
			{
				Gp = G_density(GCA,range,0);
			}
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):G",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = Gp;
			(*res)->data = (void*)result_data;
		}
			break;
		case GCALAB_ALL_PARAM:
		{
			float *result_data;
			lambdap = lambda_param(GCA);
			Zp = Z_param(GCA);
			if (samples > 0)
			{
				Gp = G_density(GCA,NULL,samples);
			}
			else
			{
				Gp = G_density(GCA,range,0);
			}
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):PA",trgt_id);
			(*res)->datalen = 3;
			result_data = (float*)malloc(3*sizeof(float));
			result_data[0] = lambdap;
			result_data[1] = Zp;
			result_data[2] = Gp;
			(*res)->data = (void*)result_data;
		}
			break;
	}
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Reverse(): compute pre-images of current CA configuration 
 */
char GCALab_OP_Reverse(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int numPreImages;
	chunk *preImages;
	GraphCellularAutomaton *GCA;
	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
	(*res) = (GCALabOutput*)malloc(sizeof(GCALabOutput)); 
	
	preImages = CAGetPreImages(GCA,&numPreImages,NULL);
	(*res)->type = CHUNK;
	sprintf((*res)->id,"(%d):R",trgt_id);
	(*res)->datalen = numPreImages*(GCA->size);
	(*res)->data = (void*)preImages;
	return GCALAB_SUCCESS;
}

/* GCALab_OP_Freq(): count CA state frequency per cell.
 */
char GCALab_OP_Freq(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int numSamples;
	unsigned int *freqs;
	int size;
	char rc;
	int i;
	chunk range[2];
	chunk *ics;
	numSamples = 0;
	for (i=0;i<nparams;i++)
	{
		if(!strcmp(params[i],"-n"))
		{
			numSamples = (int)atoi(params[++i]);
		}
		else if (!strcmp(params[i],"-l"))
		{
			range[0] = (chunk)atoi(params[++i]);
			range[1] = (chunk)atoi(params[++i]);
		}
	}
			
	GraphCellularAutomaton *GCA;
	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
	(*res) = (GCALabOutput*)malloc(sizeof(GCALabOutput)); 
				
	size = (GCA->params->N)*(GCA->params->N);
	freqs = (unsigned int *)malloc(size*sizeof(unsigned int));
	rc = GCALab_TestPointer((void*)freqs);
	if (rc <= 0)
	{
		return rc;
	}
	memset((void*)freqs,0,size*sizeof(unsigned int));

	if (numSamples)
	{
		ics = (chunk*)malloc(numSamples*sizeof(chunk));
		rc = GCALab_TestPointer((void*)ics);
		if (rc <= 0)
		{
			return rc;
		}

		for (i=0;i<numSamples;i++)
		{
			ics[i] = rand();
		}
		SumCAImages(GCA,freqs,ics,numSamples);
		free(ics);
	}
	else
	{
		SumCAImages(GCA,freqs,range,0);
	}

	(*res)->type = UINT32;
	sprintf((*res)->id,"(%d):F",trgt_id);
	(*res)->datalen = size;
	(*res)->data = (void*)freqs;
	return GCALAB_SUCCESS;
}
#else
int main(){
	printf("Install a proper operating system!\n");
	printf("See http://www.fedoraproject.org\n");
}
#endif
