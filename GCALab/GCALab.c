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
 * Last Modified: 18/01/2013
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
 *       v 0.17 (18/01/2013) - i. Extended param command to handle And average attractor
 *                                cycle length command and an average transient path
 *                                length command.
 *       v 0.18 (19/01/2013) - i. Added a command and operation registeration functions,
 *                                should have been done a while ago, but it really has 
 *                                cleaned up the GCALab Init function.
 *                             ii. removed the need to enter q-cmd to enter thread operations,
 *                                 however left it in as a command for compatability with old
 *                                 scripts.
 *                             iii. made command prompt more robust, does not seg fault on invalid
 *                                  user inputs. Meaningful error messages are returned now also.
 *                             iv. improved user input of command that require a type, no longer 
 *                                 need to know the type codes, but a name instead.
 *                             v.  imporved error handling of command so crashes are less common.
 *       v 0.19 (01/03/2012) - i. Included and tested neighbourhood type and life rule switch in
 *                                the gca create command.
 *                             ii. Added a configuration edit mode when running in Graphics mode.
 *
 * Description: Main Program for Graph Cellular Automata generation, simulation,
 *              analysis and Visualisation.
 *
 * TODO List:
 *		1. test commandline parsing - done (v 0.02)
 *		2. implement batch mode first - cancelled (v 0.04)
 *		3. use ptheads to implement a main loop plus compute threads - done (v 0.05)
 * 		4. implement text mode first, then re-implement batch mode, then graphics - done (0.12) 
 * 		5. carefully separated main engine from the mode type - done (0.11)
 * 		6. remember to comment function pointer framework carefully. 
 * 		7. should add a cycle compute function, possibly useful comparison - done (v 0.17)
 * 		8. build in more plotting/image saving options
 * 		9. add a create animation option
 * 		10. should make command more user friendly ie add a symbol table for CA and results
 * 		11. remove the need to enter q-cmd before an operation - done (0.18)
 * 		12. change input args so option id codes do not need to be known explicitly - done (0.18)
 * 		13. add register command and register operation functions - done (0.18)
 * 		14. Add cell picking so a test configuration is easy to test. Also an effective method
 * 		    of loading an initial configuration file would also be good. - done (0.19)
 * 		15. LUT colour mode need more colours (possibly a continuous colourmap) since the
 * 		    moore neighbourhood LUTs are quite large.
 * Known Issues:
 *
 *==============================================================================
 */
#define GCALAB_VERSION 0.2
#define GCALAB_AUTHOR "David J. Warne"
#define GCALAB_CW_YEAR 2015
#include "GCALab.h"
#ifdef __linux__
/*global array of workspace addresses*/
GCALab_WS **GCALab_Global;
unsigned int GCALab_numWS;
GCALab_Cmd GCALab_Cmds[GCALAB_MAXNUM_CMDS];
unsigned int GCALab_numCmds;
GCALab_Op GCALab_Ops[GCALAB_MAXNUM_OPS];
unsigned int GCALab_numOps;
unsigned char GCALab_mode;
unsigned int cur_ws;
#ifdef WITH_GRAPHICS
/* light settings*/
GLfloat ambientLight[] = { 0.1f, 0.1f, 0.1f, 1.0f };
GLfloat diffuseLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat position[] = { 10.0f, 10.0f, 10.0f, 1.0f };
GLfloat positionMirror[] = { 10.0f, -10.0f, 10.0f, 1.0f };
GLfloat fogCol[] = {0.4f,0.4f,0.4f,0.5f};
GLUquadric* quad;
float translateX,translateY,zoom;
float theta, phi;
unsigned int cur_gca;
unsigned int cur_res;
unsigned char moving;
int xprev;
int yprev;
float cellColf[8][3] = {{0.2,0.2,0.2},{1,1,1},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1}};
float lutColf[16][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1},{1,1,1},
						{1,0.5,0},{0.0,1,0.5},{0.5,0,1},{1,0.5,1},{1,1,0.5},{0.5,1,1},{0.5,0.5,0.5}};
unsigned char lutColmode;
unsigned char showMesh;
unsigned char GCALab_confedit_mode;
#endif
char stateInitials[6] = {'I','R','P','Q','E'};
char* statenames[6] = {"Idle","Running","Paused","Exiting","Error"};
char cellsymbols[8] = {' ','O','*','-','x','+','#','^'};

/**
 * @brief entry point of the GCALab Program
 * @param argc the number of commandline args
 * @param argv array of commandline args
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

/** 
 * @brief starts a GCALab session in text-only mode
 * @param opts User provided start up commandline args
 */
void GCALab_TextMode(GCALab_CL_Options* opts)
{
	char **userinput;
	int numargs;
	char rc;
	GCALab_SplashScreen();
	
	while(1)
	{
		/*get user input*/
		userinput = GCALab_CommandPrompt(&numargs);
		rc = GCALab_Process_Command(numargs,userinput);
		GCALab_HandleErr(rc);
	}
	return;
}

/**
 * @brief Runs GCALab interactively with OpenGL Grpahics
 * @param opts User provided start up commandline args
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
	showMesh = 0;
	moving = GCALAB_GRAPHICS_INTERACTION_NONE;	
	glutMainLoop();
#else
	fprintf(stderr,"Graphics Mode not enabled in this build.\n");
#endif
	return;
}

/**
 * @brief Runs Commands is batch mode, i.e., reads commands from a file
 * @param opts User provided start up commandline args
 */
char GCALab_BatchMode(GCALab_CL_Options* opts)
{
	char **userinput;
	int numargs;
	char rc;
	GCALab_SplashScreen();
	
	while(1)
	{
		/*get user input*/
		userinput = GCALab_ReadScriptCommand(opts->ScriptFile,&numargs);
		opts->ScriptFile = NULL;
		rc = GCALab_Process_Command(numargs,userinput);
		GCALab_HandleErr(rc);
	}
}

/**
 * @brief sets up GCALab with default settings, unless overridden via start-up commands
 * @details this function is also responsible for registration of extension commands
 * @param argc number of commandline args
 * @param argv array of commandline args
 * @param opts Parsed commandline options (populated by this function)
 * @returns GCALAB_SUCCESS on sucessful completion or and appropriate error code. 
 */
char GCALab_Init(int argc,char **argv,GCALab_CL_Options **opts)
{
	char *args,*desc;
	GCALab_Global = (GCALab_WS **)malloc(GCALAB_MAX_WORKSPACES*sizeof(GCALab_WS*));
	if(!(GCALab_Global))
	{
		return GCALAB_FATAL_ERROR;
	}
	
	opts[0] = GCALab_ParseCommandLineArgs(argc,argv);
	if (!(opts[0]))
	{
		return GCALAB_CL_PARSE_ERROR;
	}
	
	GCALab_numWS = 0;
	GCALab_mode = (opts[0])->mode;
	
	GCALab_numCmds = 0;
	GCALab_numOps = 0;
	
	/*register commands*/
	args = "n";
	desc = "Creates a new workspace with n CA slots.";
	GCALab_Register_Command("new-work",&GCALab_CMD_NewWorkSpace,args,desc);
	args = "none";
	desc = "Print the current workspace.";
	GCALab_Register_Command("print-work",&GCALab_CMD_PrintWorkSpace,args,desc);
	args = "none";
	desc = "Print Summary of all workspaces.";
	GCALab_Register_Command("list-work",&GCALab_CMD_ListWorkSpaces,args,desc);
	args = "id";
	desc = "Changes current workspace to be id.";
	GCALab_Register_Command("ch-work",&GCALab_CMD_ChangeWorkSpace,args,desc);
	args = "cmd";
	desc = "Enqueues the GCA operation to the current command queue.";
	GCALab_Register_Command("q-cmd",&GCALab_CMD_QueueCommand,args,desc);
	args = "cmd_id";
	desc = "Sets GCA operation to cmd_id to be ignored.";
	GCALab_Register_Command("del-cmd",&GCALab_CMD_DeleteCommand,args,desc);
	args = "none";
	desc = "The current queue will start processing.";
	GCALab_Register_Command("exec-q",&GCALab_CMD_ExecuteQueue,args,desc);
	args = "none";
	desc = "The current queue will pause after completion of the current operation.";
	GCALab_Register_Command("stop-q",&GCALab_CMD_StopQueue,args,desc);
	args = "ca_id";
	desc = "Prints GCA with id ca_id in the current workspace.";
	GCALab_Register_Command("print-ca",&GCALab_CMD_PrintCA,args,desc);
	args = "ca_id";
	desc = "Prints GCA evolution.";
	GCALab_Register_Command("print-st",&GCALab_CMD_PrintSTP,args,desc);
	args = "res_id";
	desc = "Prints result data with id res_id.";
	GCALab_Register_Command("print-res",&GCALab_CMD_PrintResults,args,desc);
	args = "none";
	desc = "Exit GCALab.";
	GCALab_Register_Command("quit",&GCALab_CMD_Quit,args,desc);
	args = "none";
	desc = "Prints this help menu.";
	GCALab_Register_Command("help",&GCALab_CMD_PrintHelp,args,desc);
	args = "none";
	desc = "Prints available compute commands.";
	GCALab_Register_Command("list-cmds",&GCALab_CMD_PrintOperations,args,desc);
	/*register operations*/	
	args = "none";
	desc = "No Operation";
	GCALab_Register_Operation("nop",&GCALab_OP_NOP,args,desc);
	args = "i -f filename";
	desc = "Loads a *.gca file into the current workspace";
	GCALab_Register_Operation("load",&GCALab_OP_Load,args,desc);
	args = "i -f filename (-g | -r)";
	desc = "Save a GCA or result with id to file";
	GCALab_Register_Operation("save",&GCALab_OP_Save,args,desc);
	args = "i -t Tfinal [-I] [-f icfile | -c (random | point | checker | stripe)]";
	desc = "simulates the id to Tfinal";
	GCALab_Register_Operation("sim",&GCALab_OP_Simulate,args,desc);
	args = "i (((-m meshfile | -t numcells genus) -s numstates -r (code | totalistic | thresh | life ) rulecode) | -eca numCells numNeighbours rulecode) [-c (random | point | checker | stripe)] [-w windowsize] [-nh (neumann | moore)]";
	desc = "Creates a new graph cellular automaton in the current workspace";
	GCALab_Register_Operation("gca",&GCALab_OP_GCA,args,desc);
    args = "i [-p prob]";
    desc = "Rotate neighbourhoods with probability p";
	GCALab_Register_Operation("rotate",&GCALab_OP_Rotate,args,desc);
	args = "i -n numsamples -t timesteps -e entropytype -p";
	desc = "Computes entropy measures of graph cellular automaton at i";
	GCALab_Register_Operation("entropy",&GCALab_OP_Entropy,args,desc);
	args = "i -p paramtype [-l config0 configN | -n numSamples -t maxT]";
	desc = "Computes complexity parameters such as Langton's lambda";
	GCALab_Register_Operation("param",&GCALab_OP_Param,args,desc);
	args = "i";
	desc = "Computes pre-images of the current configuration of the graph cellular automaton at i";
	GCALab_Register_Operation("pre",&GCALab_OP_Reverse,args,desc);
	args = "i (-n numsamples | -l config0 configN)";
	desc = "Computes state frequency histogram for each cell in the graph cellular automaton at i";
	GCALab_Register_Operation("freq",&GCALab_OP_Freq,args,desc);
	args = "i -t timesteps";
	desc = "Computes the non-quiescient population density over time.";
	GCALab_Register_Operation("pop",&GCALab_OP_Pop,args,desc);
	return GCALAB_SUCCESS;
}

/**
 * @brief Registers command functions for the GCALab command prompt.
 * @param id the name of the command.
 * @param f function pointer to handle the command.
 * @param args a human readable list of arguments.
 * @param desc a human readable summary of the function.
 */
void GCALab_Register_Command(char *id,char (*f)(int,char**),char * args, char * desc)
{
	if (GCALab_numCmds < GCALAB_MAXNUM_CMDS)
	{
		GCALab_Cmds[GCALab_numCmds].id = id;
		GCALab_Cmds[GCALab_numCmds].f = f;
		GCALab_Cmds[GCALab_numCmds].args = args;
		GCALab_Cmds[GCALab_numCmds].desc = desc;
		GCALab_numCmds++;
	}
}

/** 
 * @brief Registers a compute operation
 * @param id the name of the operation.
 * @param f function pointer to handle the operation.
 * @param args a human readable list of arguments.
 * @param desc a human readable summary of the operation.
 */
void GCALab_Register_Operation(char *id,char (*f)(unsigned char,unsigned int,int,char**,GCALabOutput**),char* args,char * desc)
{
	if (GCALab_numOps < GCALAB_MAXNUM_OPS)
	{
		GCALab_Ops[GCALab_numOps].id = id;
		GCALab_Ops[GCALab_numOps].f = f;
		GCALab_Ops[GCALab_numOps].args = args;
		GCALab_Ops[GCALab_numOps].desc = desc;
		GCALab_numOps++;
	}
}

/**
 * @brief Interprets the user command prompt input as returned by GCALab_CommandPrompt()
 * @param nargs the number of arguments
 * @param args the argument list
 * @returns the return code of the user selected function
 */
char GCALab_Process_Command(int nargs,char **args)
{
	char rc;
	int i;
	char *cmd;
	/*ensure input is valid*/
	rc = GCALab_TestPointer((void*)args);
	GCALab_HandleErr(rc);
		
	cmd = args[0];
	for (i=0;i<GCALab_numCmds;i++)
	{
		if (!strcmp(cmd,GCALab_Cmds[i].id))
		{
			rc = (*(GCALab_Cmds[i].f))(nargs,args);
			GCALab_HandleErr(rc);
			break;
		}
	}
	/*if we don't find a command then assume a q-cmd*/
	/*and pass it down to the queue*/
	if (i == GCALab_numCmds)
	{
		rc = GCALab_CMD_QueueCommand(nargs,args);	
	}

	/*clean up user args*/
	for (i=0;i<nargs;i++)
	{
		free(args[i]);
	}
	free(args);

	return rc;
}

/**
 * @brief returns an error if the pointer is null
 * @param ptr the pointer to test.
 * @retVal GCALAB_SUCCESS if ptr != NULL
 * @retVal GCALAB_MEM_ERROR if ptr == NULL
 */
char GCALab_TestPointer(void* ptr)
{
	if (!ptr)
		return GCALAB_MEM_ERROR;
	else
		return GCALAB_SUCCESS;
}

/**
 * @brief tests if the given workspace id maps to a real workspace
 * @param ws_id the current workspace id.
 */
char GCALab_ValidWSId(unsigned int ws_id)
{
	/*for now simple*/
	if (ws_id >= GCALab_numWS)
		return GCALAB_INVALID_WS_ERROR;
	else
		return GCALAB_SUCCESS;
}

/**
 * @brief tests if the given returen code indicates an error,if the error is fatal then it 
 * aborts.
 * @param rc the return code to test
 */
void GCALab_HandleErr(char rc)
{
	if (rc <= 0)
	{
		switch (rc)
		{
			case GCALAB_MEM_ERROR:
				fprintf(stderr,"OUT OF MEMORY!!!\n");
				GCALab_HandleErr(GCALAB_FATAL_ERROR);
				break;
			case GCALAB_FATAL_ERROR:
				fprintf(stderr,"A fatal error occurred. Aborting.\n");
				exit(1);
				break;
			case GCALAB_INVALID_OPTION:
				fprintf(stderr,"Invalid command option.\n");
				break;
			case GCALAB_UNKNOWN_OPTION:
				fprintf(stderr,"Unknown command.\n");
				break;
			case GCALAB_CL_PARSE_ERROR:
				break;
			case GCALAB_INVALID_WS_ERROR:
				fprintf(stderr,"Invalid Workspace Id.\n");
				break;
			case GCALAB_THREAD_ERROR:
				fprintf(stderr,"Something funky happened with a thread.\n");
				GCALab_HandleErr(GCALAB_FATAL_ERROR);
				break;
		}
	}
}

/**
 * @brief converts the command string into a command id
 * @param cmd command name to lookup the op-code of.
 * @returns the command op-code or GCALAB_NOP if not found.
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

/**
 * @brief Processing thread attached to a workspace.
 * @param params thread args, just the workspace id.
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
		/*Worker thread finite state machine*/
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

/**
 * @brief evaluates to true when the Command queue for the given workspace is empty.
 * @param ws_id the workspace id to test.
 * @returns 1 if the commandqueue is empty, 0 otherwise.
 */
unsigned char GCALab_CommandQueueEmpty(unsigned char ws_id)
{
	unsigned char empty;
	GCALab_LockWS(ws_id);
	empty = (WS(ws_id)->numcommands == 0);
	GCALab_UnLockWS(ws_id);
	return empty;
}

/**
 * @brief executes the next command in the given command queue.
 * @param ws_id the workspace id to process
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
/**
 * @brief Prints Author and affiliation information
 */
void GCALab_PrintAbout(void)
{
	fprintf(stdout,"\n The Graph Cellular Automata Lab is a multi-threaded,\n");
	fprintf(stdout," flexible, analysis tool forthe exploration of cellular automata\n"); 
	fprintf(stdout," defined on graphs.\n");
	fprintf(stdout,"\n Version:\t%0.2f\n",GCALAB_VERSION);
	fprintf(stdout," Author:\t%s\n",GCALAB_AUTHOR);
	fprintf(stdout," School:\tSchool of Electrical Engineering and Computer Science,\n");
	fprintf(stdout," \t\tThe Queensland University of Technology, Brisbane, Australia\n");
	fprintf(stdout," Contact:\tdavid.warne@qut.edu.au\n");
	fprintf(stdout," url:\thttps://github.com/davidwarne/GCALab.git\n");
	return;
}

/**
 * @brief prints user license summary
 */
void GCALab_PrintLicense(void)
{
	fprintf(stdout,"\n GCALab v %0.2f  Copyright (C) %d  %s\n",GCALAB_VERSION,GCALAB_CW_YEAR,GCALAB_AUTHOR);
	fprintf(stdout," This program comes with ABSOLUTELY NO WARRANTY.\n");
	fprintf(stdout," This is free software, and you are welcome to redistribute it\n");
	fprintf(stdout," under certain conditions.\n");
	return;
}

/**
 * @brief startup screen
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
			GCALab_PrintLicense();
			break;
	}
	return;
}

/**
 * @brief Prints the help menu
 */
void GCALab_PrintUsage()
{
	printf("CGALab - Options");
	printf("\t [-g,--graphics]\n\t\t : start in interactive mode\n");
	printf("\t [-i,--interactive]\n\t\t : start in interactive mode\n");
	printf("\t [-b,--batch]\n\t\t : start in batch mode\n");
}

/**
 * @brief Prompts the user to enter a command in text mode, 
 * @details the users input string is then tokenised with whitespace as the delimiter.
 * @param numargs populated with the numerb of arguments read
 * @returns the array of arguments
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

/**
 * @brief Prompts the user to enter a command in text mode, 
 * @details the users input string is then tokenised with whitespace as the delimiter.
 * @param filename the name of the script file
 * @param numargs populated with the numerb of arguments read
 * @returns the array of arguments
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

/**
 * @brief Initialises the options struct with defaults values
 * @param opts user commandline options
 */
void GCALab_InitCL_Options(GCALab_CL_Options* opts)
{
	opts->mode = GCALAB_DEFAULT_MODE;
	opts->load_CA = 0;
	opts->save_CA = 0;
	opts->CAInputFilename = "";
	opts->CAOutputFilename = "";
}

/**
 * @brief parses the command line arguments and returns a struct of user options.
 * @param argc number of commandline args as proided from main
 * @param argv the array of commandline args 
 * @returns a options structure with fields populated as per the input args.
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
			else /*unknown option*/
			{
				GCALab_PrintUsage();
				return NULL;
			}
			
		}
	}
	
	return CL_opt;	
}

/**
 * @brief Creates a new processing workspace, 
 * @details the user can set the limit on the number of CA objects the workspace can hold.
 * @param GCALimit the maximum number of CA that can be handled by this workspace
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

/**
 * @brief aquire a lock on the given workspace
 * @param ws_id the workspace to lock
 */
void GCALab_LockWS(unsigned char ws_id)
{
	pthread_mutex_lock(&(GCALab_Global[ws_id]->wslock));
}

/**
 * @brief release a lock on the given workspace
 * @parakm ws_id the workspace to unlock
 */
void GCALab_UnLockWS(unsigned char ws_id)
{
	pthread_mutex_unlock(&(GCALab_Global[ws_id]->wslock));
}

/**
 * @brief appends the given command string to the command queue of the given workspace.
 * @param ws_id workspace to modify
 * @param command_id the op-code of the user queued command
 * @param target_id the location of the results in the data array
 * @param params list of args to the user command
 * @param nparams the number of args to teh user command
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

/**
 * @brief cancels the the command located at address index of the command 
 * queue of the given workspace.
 * @param ws_id the workspace we are modifying
 * @param index the index of the command in the queue
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

/**
 * @brief Tells workspace to continue processing
 * @param ws_id the workspace to modify
 */
char GCALab_ProcessCommandQueue(unsigned char ws_id)
{
	GCALab_LockWS(ws_id);
	WS(ws_id)->state = GCALAB_WS_STATE_IDLE;
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/**
 * @brief Tells workspace to halt processing
 * @param ws_id the workspace to modify
 */
char GCALab_PauseCommandQueue(unsigned char ws_id)
{
	GCALab_LockWS(ws_id);
	WS(ws_id)->state = GCALAB_WS_STATE_PAUSED;
	GCALab_UnLockWS(ws_id);
	return GCALAB_SUCCESS;
}

/**
 * @brief Tells workspace print all queued commands
 * @param ws_id the workspace to modify
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

/**
 * @brief gets the state flag of a workspace
 * @param ws_id the workspace to modify
 * @returns the state of the workspace
 */
unsigned int GCALab_GetState(unsigned char ws_id)
{
	unsigned int state;
	GCALab_LockWS(ws_id);
	state = GCALab_Global[ws_id]->state;
	GCALab_UnLockWS(ws_id);
 	return state;
}

/**
 * @brief sets the state flag of a workspace safely
 * @param ws_id the workspace to modify
 * @param state the new state of the workspace
 */
void GCALab_SetState(unsigned char ws_id,unsigned int state)
{
	GCALab_LockWS(ws_id);
	GCALab_Global[ws_id]->state = state;
	GCALab_UnLockWS(ws_id);
	return;
}

/**
 * @brief Nicely and humainly kills the session
 * @param rc return code
 * @note no real cellular automata where wounded in this process :)
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
			usleep(1);
		}

	}
	for (i=0;i<GCALab_numWS;i++)
	{
		GCALab_SetState(i,GCALAB_WS_STATE_EXITING);
	}
	exit(0);
}

/**
 * @brief Prints the help menu
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

/**
 * @brief Print list of queuable compute functions
 */
void GCALab_PrintOperations(void)
{
	int i;

	fprintf(stdout,"\nOperations List:\n");
	fprintf(stdout,"-----------------\n");
	for (i=0;i<GCALab_numOps;i++)
	{
		fprintf(stdout,"\tSynopsis: %s %s\n",GCALab_Ops[i].id,GCALab_Ops[i].args);
		fprintf(stdout,"\tDescription: %s\n\n",GCALab_Ops[i].desc);
	}
	return;
}

/**
 * @brief prints summart information about the given workspace.
 * @param ws_id the workspace to print
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
	glLightfv(GL_LIGHT0, GL_POSITION, position); 
	/*smooth shading model*/
	glShadeModel(GL_FLAT);
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

	glPushMatrix();
	/* draw the reflected scene first*/
	glLightfv(GL_LIGHT0,GL_POSITION,positionMirror);
	/*transformations*/
	glTranslatef(translateX,translateY,zoom-10.5f);
	glRotatef(phi+30.0f,1.0f,0.0f,0.0f);
	glRotatef(theta,0.0f,1.0f,0.0f);
	
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
    glPopMatrix();

	glPushMatrix();
	/* draw the "real" scene... but how do we define real?*/
	glLightfv(GL_LIGHT0, GL_POSITION, position); 
	/*transformations*/
	glTranslatef(translateX,translateY,zoom-10.5f);
	glRotatef(phi+30.0f,1.0f,0.0f,0.0f);
	glRotatef(theta,0.0f,1.0f,0.0f);
	GCALab_Graphics_DrawGrid();

	
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
	glPopMatrix();
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
	
	if(!(n = (float*)malloc(3*sizeof(float))))
	{
		return;
	}

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
		n = Normal_f(m->vList->verts + 3*(m->fList->faces[3*i]),m->vList->verts + 3*(m->fList->faces[3*i+1]),m->vList->verts + 3*(m->fList->faces[3*i+2]),n);
		/*draw the cell*/
		glNormal3fv(n);
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i]));
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+1]));
		glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+2]));
	}
	glEnd();
	
    free(n);
	if (showMesh)
	{
    	glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glColor4f(0.0,0.0,0.0,1.0);
		for (i=0;i<m->fList->numFaces;i++)
		{
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i]));
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+1]));
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+1]));
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+2]));
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+2]));
			glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i]));
		}
		glEnd();
    	glEnable(GL_LIGHTING);
	}
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

/* GCALab_Graphics_PickCell(): to be used with object picking in 
 * configuration ediut mode
 */
unsigned int GCALab_Graphics_PickCell(int x,int y)
{
	mesh *m;
	unsigned int i;
	GraphCellularAutomaton *gca;
	float f;
	GLint v[4];
	GLubyte tag[3];
	glGetIntegerv(GL_VIEWPORT,v);
	/*clear screen for a new render*/
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	
	/*transformations*/
	glTranslatef(translateX,translateY,zoom-10.5f);
	glRotatef(phi+30.0f,1.0f,0.0f,0.0f);
	glRotatef(theta,0.0f,1.0f,0.0f);
	

	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_MULTISAMPLE); 
	glDisable(GL_FOG);
	glPushMatrix();
	glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
	glPopMatrix();
	if (GCALab_numWS > 0)
	{
		if (WS(cur_ws)->numGCA > 0)
		{
			glPushMatrix();
			glTranslatef(0.0,GCALAB_GRAPHICS_OFFSET,0.0);
			m = WS(cur_ws)->GCAGeometry[cur_gca];
			gca = WS(cur_ws)->GCAList[cur_gca];
	
			glBegin(GL_TRIANGLES);
			for (i=0;i<m->fList->numFaces;i++)
			{
				tag[2] = i & 0xFF;
				tag[1] = (i >> 8) & 0xFF;
				tag[0] = (i >> 16) & 0xFF;
				/*change material properties based on cell state*/
				glColor3ub(tag[0],tag[1],tag[2]);
				glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i]));
				glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+1]));
				glVertex3fv(m->vList->verts + 3*(m->fList->faces[3*i+2]));
			}
			glEnd();
			glPopMatrix();
		}
	} 
	glReadPixels(x,v[3]-1-y,1,1,GL_RGB,GL_UNSIGNED_BYTE,(void*)tag);
	i = 0;
	glEnable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_FOG);
	i = (((unsigned int)tag[0]) << 16) | (((unsigned int)tag[1]) << 8) | ((unsigned int)tag[2]);
	return i;
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
			/*get user input*/
			userinput = GCALab_CommandPrompt(&numargs);
			rc = GCALab_Process_Command(numargs,userinput);
			GCALab_HandleErr(rc);
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
		case 'p':
			GCALab_confedit_mode = !GCALab_confedit_mode;
			break;
		case 'm':
			showMesh = !showMesh;
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
	if (!GCALab_confedit_mode)
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
	else /*configuration edit mode*/
	{
		/*only continue if the id was valid*/
		if (GCALab_ValidWSId(cur_ws) != GCALAB_INVALID_WS_ERROR)
		{
			if (WS(cur_ws)->numGCA != 0)
			{
				switch(button)
				{
					case GLUT_LEFT_BUTTON:
						SetCellStatePacked(WS(cur_ws)->GCAList[cur_gca],GCALab_Graphics_PickCell(x,y),1);
						break;
					case GLUT_RIGHT_BUTTON:
						SetCellStatePacked(WS(cur_ws)->GCAList[cur_gca],GCALab_Graphics_PickCell(x,y),0);
						break;
				}
			}
		}
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
	char rc;
	rc = GCALab_ValidWSId(cur_ws);
	/*only continue if the id was valid*/
	if (rc == GCALAB_INVALID_WS_ERROR) return rc;
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
		int i,s;
		/*first get the id map for the command*/
		s = 1;
		cmd_code = GCALab_GetCommandCode(argv[s]);
		if (cmd_code == GCALAB_NOP)
		{
			s = 0;
			/*try the first arg (incase q-cmd is not included)*/
			cmd_code = GCALab_GetCommandCode(argv[s]);
			if (cmd_code == GCALAB_NOP)
			{
				return GCALAB_UNKNOWN_OPTION;
			}
		}
		/*the target id is next*/
		target = (unsigned int)atoi(argv[s+1]);
			
		/*everything else is specific to the command*/
		numparams = argc - (s+2);
		/*make a copy*/
		params = strvncpy(argv+(s+2),numparams,GCALAB_MAX_STRLEN);
		rc = GCALab_TestPointer((void*)params);
		if (rc == GCALAB_MEM_ERROR) return rc;
		rc = GCALab_ValidWSId(cur_ws);
		if (rc < 0) return rc;
		/*push it to the queue*/
		return GCALab_QueueCommand(cur_ws,cmd_code,target,params,numparams);	
	}
}

/* GCALab_CMD_DeleteCommand(): GCALab command to delete a compute
 *                             operation.
 */
char GCALab_CMD_DeleteCommand(int argc, char **argv)
{
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int ind;
		ind = (unsigned int)atoi(argv[1]);
		rc = GCALab_ValidWSId(cur_ws);
		/*only continue if the id was valid*/
		if (rc == GCALAB_INVALID_WS_ERROR) return rc;
		return GCALab_CancelCommand(cur_ws,ind);
	}
}

/* GCALab_CMD_ExecuteQueue(): GCALab Command to signal a thread to 
 *                            process a command queue
 */
char GCALab_CMD_ExecuteQueue(int argc, char **argv)
{
	char rc;
	rc = GCALab_ValidWSId(cur_ws);
	/*only continue if the id was valid*/
	if (rc == GCALAB_INVALID_WS_ERROR) return rc;
	return GCALab_ProcessCommandQueue(cur_ws);
}

/* GCALab_CMD_StopQueue(): GCALab command to signal a thread to
 *                         halt processing of a command queue.
 */
char GCALab_CMD_StopQueue(int argc, char **argv)
{
	char rc;
	rc = GCALab_ValidWSId(cur_ws);
	/*only continue if the id was valid*/
	if (rc == GCALAB_INVALID_WS_ERROR) return rc;
	return GCALab_PauseCommandQueue(cur_ws);
}

/* GCALab_CMD_PrintCA(): GCALab command to print CA summary
 */
char GCALab_CMD_PrintCA(int argc, char **argv)
{
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		rc = GCALab_ValidWSId(cur_ws);
		/*only continue if the id was valid*/
		if (rc == GCALAB_INVALID_WS_ERROR) return rc;
		return GCALab_PrintCA(cur_ws,id);
	}
}

/* GCALab_CMD_PrintSTP(): GCALab command to print CA evolution
 */
char GCALab_CMD_PrintSTP(int argc, char ** argv)
{
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		rc = GCALab_ValidWSId(cur_ws);
		/*only continue if the id was valid*/
		if (rc == GCALAB_INVALID_WS_ERROR) return rc;
		return GCALab_PrintSTP(cur_ws,id);
	}
}

/* GCALab_CMD_PrintResults(); GCALab command to print result data.
 */
char GCALab_CMD_PrintResults(int argc, char **argv)
{
	char rc;
	if (argc < 2)
	{
		return GCALAB_INVALID_OPTION;
	}
	else
	{
		unsigned int id;
		id = (unsigned int)atoi(argv[1]);
		rc = GCALab_ValidWSId(cur_ws);
		/*only continue if the id was valid*/
		if (rc == GCALAB_INVALID_WS_ERROR) return rc;
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
			char * typestr = params[++i];
			if (!strcmp(typestr,"point"))
			{
				ic_type = POINT_IC_TYPE;
			}
			else if (!strcmp(typestr,"random"))
			{
				ic_type = NOISE_IC_TYPE;
			}
			else if (!strcmp(typestr,"stripe"))
			{
				ic_type = STRIPE_IC_TYPE;
			}
			else if (!strcmp(typestr,"checker"))
			{
				ic_type = CHECKER_IC_TYPE;
			}
			else
			{
				return GCALAB_INVALID_OPTION;
			}
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
	unsigned char r_type,nh_type;
	unsigned char s,k,eca,ic_type;
	int  i;
	char rc;
	GraphCellularAutomaton *GCA;
	CellularAutomatonParameters *CAparams;
	mesh *m;

	windowsize = 0;
	meshfile = NULL;
	eca = 0;
	nh_type = DEFAULT_NEIGHBOURHOOD_TYPE;
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
			char * typestr = params[++i];
			if (!strcmp(typestr,"code"))
			{
				r_type = CODE_RULE_TYPE;
			}
			else if (!strcmp(typestr,"thresh"))
			{
				r_type = THRESH_RULE_TYPE;
			}
			else if (!strcmp(typestr,"totalistic"))
			{
				r_type = COUNT_RULE_TYPE;
			}
			else if (!strcmp(typestr,"life"))
			{
				r_type = LIFE_RULE_TYPE;
			}
			else
			{
				return GCALAB_INVALID_OPTION;
			}
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
			char * typestr = params[++i];
			if (!strcmp(typestr,"point"))
			{
				ic_type = POINT_IC_TYPE;
			}
			else if (!strcmp(typestr,"random"))
			{
				ic_type = NOISE_IC_TYPE;
			}
			else if (!strcmp(typestr,"stripe"))
			{
				ic_type = STRIPE_IC_TYPE;
			}
			else if (!strcmp(typestr,"checker"))
			{
				ic_type = CHECKER_IC_TYPE;
			}
			else
			{
				return GCALAB_INVALID_OPTION;
			}
		}
		else if (!(strcmp(params[i],"-nh")))
		{
			char * typestr = params[++i];
			if (!strcmp(typestr,"neumann"))
			{
				nh_type = VON_NEUMANN_NEIGHBOURHOOD_TYPE;
			}
			else if (!strcmp(typestr,"moore"))
			{
				nh_type = MOORE_NEIGHBOURHOOD_TYPE;
			}
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
		CAparams = CreateCAParams(nh_type,m,s,r_type,r,windowsize);
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

/* GCALab_OP_Rotate(): rotate cell neighbourhoods at random
 */
char GCALab_OP_Rotate(unsigned char ws_id,unsigned int trgt_id, int nparams, char ** params, GCALabOutput **res)
{
    int i;
    float p;
    GraphCellularAutomaton *GCA;
    p = 0.5;
	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
    for (i=0;i<nparams;i++)
    {
        if (!strcmp(params[i],"-p"))
        {
            p = (float)atof(params[++i]);
        }
        else 
        {
            return GCALAB_INVALID_OPTION;
        }
    }

    for (i=0;i<GCA->params->N;i++)
    {
        if (((float)rand())/((float)RAND_MAX) > p)
        {
            RotateNeighbourhood(GCA,i,rand());
        }
    }
	
    return GCALAB_SUCCESS;
}
/* GCALab_OP_Entropy(): compute entropy meansures on a given CA
 */
char GCALab_OP_Entropy(unsigned char ws_id,unsigned int trgt_id,int nparams, char ** params,GCALabOutput **res)
{
	unsigned int numSamples;
	unsigned int T;
	unsigned type;
    unsigned int rotate;
	int i,j;
	char rc;
	GraphCellularAutomaton *GCA;
	
	float *p,*logs_p, *S_i,*pt,*logs_pt,*IE,*logQ;
	unsigned char *flags;
	unsigned int *count,*countt,*wl,*Q;
	
	float S_mu,W_mu,I_mu,I_sigma;
	float *result_data;
	numSamples = 1;
    rotate = 0;
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
        else if (!strcmp(params[i],"-p"))
        {
            rotate = 1;
        }
		else if(!strcmp(params[i],"-e"))
		{
			char * typestr = params[++i];
			if (!strcmp(typestr,"Shannon"))
			{
				type = GCALAB_SHANNON_ENTROPY;
			}
			else if (!strcmp(typestr,"Word"))
			{
				type = GCALAB_WORD_ENTROPY;
			}
			else if (!strcmp(typestr,"Input"))
			{
				type = GCALAB_INPUT_ENTROPY;
			}
			else if (!strcmp(typestr,"All"))
			{
				type = GCALAB_ALL_ENTROPY;
			}
			else
			{
				return GCALAB_INVALID_OPTION;
			}
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
            if (rotate)
            {
                for (i=0;i<numSamples;i++)
                {
                    for (j=0;j<GCA->params->N;j++)
                    {
                        RotateNeighbourhood(GCA,j,rand());
                    }
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			    	S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
            
                }
            }
            else
            {
		        for (i=0;i<numSamples;i++)
			    {
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			    	S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
			    }
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
            if (rotate)
            {
			    for (i=0;i<numSamples;i++)
			    {
                    for (j=0;j<GCA->params->N;j++)
                    {
                        RotateNeighbourhood(GCA,j,rand());
                    }
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			        W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			    }
            }
            else
            {
			    for (i=0;i<numSamples;i++)
			    {
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			        W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			    }
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
            if (rotate)
            {
                for (j=0;j<GCA->params->N;j++)
                {
                    RotateNeighbourhood(GCA,j,rand());
                }
			}
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
            if (rotate)
            {
                for (i=0;i<numSamples;i++)
                {
                    for (j=0;j<GCA->params->N;j++)
                    {
                        RotateNeighbourhood(GCA,j,rand());
                    }
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			    	S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
            
                }
            }
            else
            {
		        for (i=0;i<numSamples;i++)
			    {
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			    	S_mu += ShannonEntropy(GCA,T,p,logs_p,S_i,flags,count);
			    }
            }
			S_mu = S_mu/((float)numSamples);
					
			/*compute avg word entropy*/
            if (rotate)
            {
			    for (i=0;i<numSamples;i++)
			    {
                    for (j=0;j<GCA->params->N;j++)
                    {
                        RotateNeighbourhood(GCA,j,rand());
                    }
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			        W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			    }
            }
            else
            {
			    for (i=0;i<numSamples;i++)
			    {
			    	ResetCA(GCA);
			    	SetCAIC(GCA,NULL,NOISE_IC_TYPE);
			        W_mu += WordEntropy(GCA,T,pt,logs_pt,S_i,flags,countt,wl);
			    }
            }
			W_mu = W_mu/((float)numSamples);
					
			/*compute avg and varience of I*/
            if (rotate)
            {
                for (j=0;j<GCA->params->N;j++)
                {
                    RotateNeighbourhood(GCA,j,rand());
                }
			}
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
	float lambdap,Zp,Gp,Cp,Tp;
	int i;
	unsigned int samples,maxT;
	GraphCellularAutomaton *GCA;
	samples = 0;
	maxT = 1200;
	for (i=0;i<nparams;i++)
	{
		if (!strcmp(params[i],"-p"))
		{
			char *typestr = params[++i];
			if (!strcmp(typestr,"lambda"))
			{
				type = GCALAB_LAMBDA_PARAM;
			}
			else if (!strcmp(typestr,"Z"))
			{
				type = GCALAB_Z_PARAM;
			}
			else if (!strcmp(typestr,"G-density"))
			{
				type = GCALAB_G_PARAM;
			}
			else if (!strcmp(typestr,"Att-length"))
			{
				type = GCALAB_C_PARAM;
			}
			else if (!strcmp(typestr,"Trans-length"))
			{
				type = GCALAB_T_PARAM;
			}
			else if (!strcmp(typestr,"All"))
			{
				type = GCALAB_ALL_PARAM;
			}
			else
			{
				return GCALAB_INVALID_OPTION;
			}
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
		else if (!strcmp(params[i],"-t"))
		{
			maxT = (unsigned int)atoi(params[++i]);
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
			/*compute Langton's lambda parameter*/
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
			/*compute Wuensche's Z parameter*/
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
			/*compute the density of Garden-of-Eden configurations*/
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
		case GCALAB_C_PARAM:
		{
			float *result_data;
			/*compute the average attractor cycle length*/
			if (samples > 0)
			{
				Cp = AttLength(GCA,NULL,samples,maxT);
			}
			else
			{
				Cp = AttLength(GCA,range,0,maxT);
			}
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):C",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = Cp;
			(*res)->data = (void*)result_data;
		}
			break;
		case GCALAB_T_PARAM:
		{
			float *result_data;
			/*compute the average transient path length*/
			if (samples > 0 )
			{
				Tp = TransLength(GCA,NULL,samples,maxT);
			}
			else
			{
				Tp = TransLength(GCA,range,0,maxT);
			}
			/*store outputs*/
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):T",trgt_id);
			(*res)->datalen = 1;
			result_data = (float*)malloc(sizeof(float));
			result_data[0] = Tp;
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
				Cp = AttLength(GCA,NULL,samples,maxT);
				Tp = TransLength(GCA,NULL,samples,maxT);
			}
			else
			{
				Gp = G_density(GCA,range,0);
				Cp = AttLength(GCA,range,0,maxT);
				Tp = TransLength(GCA,range,0,maxT);
			}
			(*res)->type = FLOAT32;
			sprintf((*res)->id,"(%d):PA",trgt_id);
			(*res)->datalen = 5;
			result_data = (float*)malloc(((*res)->datalen)*sizeof(float));
			result_data[0] = lambdap;
			result_data[1] = Zp;
			result_data[2] = Gp;
			result_data[3] = Cp;
			result_data[4] = Tp;
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
	GraphCellularAutomaton *GCA;
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

char GCALab_OP_Pop(unsigned char ws_id,unsigned int trgt_id,int argc, char ** argv,GCALabOutput **res)
{
	unsigned int i,T;
	char rc;
	float * dense;
	GraphCellularAutomaton *GCA;
	
	for (i=0;i<argc;i++)
	{
		if (!strcmp(argv[i],"-t"))
		{
			T = (unsigned int)atoi(argv[++i]);
		}
	}

	/*Grab a reference to the CA we want to play with*/
	GCA = WS(ws_id)->GCAList[trgt_id];
	(*res) = (GCALabOutput*)malloc(sizeof(GCALabOutput)); 
				
	dense = (float *)malloc(T*sizeof(float));
	rc = GCALab_TestPointer((void*)dense);
	if (rc <= 0)
	{
		return rc;
	}

	PopDensity(GCA,NULL,T,dense);

	(*res)->type = FLOAT32;
	sprintf((*res)->id,"(%d):D",trgt_id);
	(*res)->datalen = T;
	(*res)->data = (void*)dense;
	
	return GCALAB_SUCCESS;
}
#else
int main(){
	printf("Install a proper operating system!\n");
	printf("See http://www.fedoraproject.org\n");
}
#endif
