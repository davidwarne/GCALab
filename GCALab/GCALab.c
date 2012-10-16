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
 * Known Issues:
 *
 *==============================================================================
 */

#include "GCALab.h"

/*global array of workspace addresses*/
GCALab_WS **GCALab_Global;
unsigned int GCALab_numWS;
unsigned char GCALab_mode;


/* main(): entry point of the GCALab Program
 */
int main(int argc, char **argv)
{
	GCALab_CL_Options* CL_opt;
	char rc;
	
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
	unsigned int cur_ws;
	GCALab_SplashScreen();
	
	while(1)
	{
		/*get user input*/
		userinput = GCALab_CommandPrompt(&numargs);
		/*ensure input is valid*/
		rc = GCALab_TestPointer((void*)userinput);
		GCALab_HandleErr(rc);
		
		cmd = userinput[0];

		if (!strcmp(cmd,"nwork"))
		{
			int lim;
			if (numargs < 2)
			{
				GCALab_HandleErr(GCALAB_INVALID_OPTION);
			}
			else
			{
				lim = atoi(userinput[1]);
				rc = GCALab_NewWorkSpace(lim);
				GCALab_HandleErr(rc);
				cur_ws = GCALab_numWS - 1;
				printf("New Workspace create! ID = %d\n",cur_ws);
			}
		}
		else if(!strcmp(cmd,"printwork"))
		{
			rc = GCALab_PrintWorkSpace(cur_ws);
			GCALab_HandleErr(rc);
		}
		else if (!strcmp(cmd,"listallwork"))
		{
			rc = GCALab_ListWorkSpaces();
			GCALab_HandleErr(rc);
		}
		else if (!strcmp(cmd,"help"))
		{
			GCALab_PrintHelp();
		}
		else if (!strcmp(cmd,"chwork"))
		{
			if (numargs < 2)
			{
				GCALab_HandleErr(GCALAB_INVALID_OPTION);
			}
			else
			{
				unsigned int new_ws;
				new_ws = (unsigned int)atoi(userinput[1]);
				rc = GCALab_ValidWSId(new_ws);
				GCALab_HandleErr(rc);
				/*only update if the id was valid*/
				if (rc != GCALAB_INVALID_WS_ERROR)
				{
					cur_ws = new_ws;
				}
				fprintf(stdout,"Current Workspace is ID = %d\n",cur_ws);
			}
		}
		else if (!strcmp(cmd,"q"))
		{
			if (numargs < 3)
			{
				GCALab_HandleErr(GCALAB_INVALID_OPTION);
			}
			else
			{
				unsigned int cmd_code;
				unsigned int target;
				int numparams;
				char ** params;
				int i;
				/*first get the id map for the command*/
				cmd_code = GCALab_GetCommandCode(userinput[1]);
				/*the target id is next*/
				target = (unsigned int)atoi(userinput[2]);
			
				/*everything else is specific to the command*/
				numparams = numargs - 3;
				/*make a copy*/
				params = strvncpy(userinput+3,numparams,GCALAB_MAX_STRLEN);
				rc = GCALab_TestPointer((void*)params);
				GCALab_HandleErr(rc);
			
				/*push it to the queue*/
				rc = GCALab_QueueCommand(cur_ws,cmd_code,target,params,numparams);	
				GCALab_HandleErr(rc);
			}
		}
		else if (!strcmp(cmd,"cancel"))
		{
			
			if (numargs < 2)
			{
				GCALab_HandleErr(GCALAB_INVALID_OPTION);
			}
			else
			{
				unsigned int ind;
				ind = (unsigned int)atoi(userinput[1]);
				rc = GCALab_CancelCommand(cur_ws,ind);
				GCALab_HandleErr(rc);
			}
		}
		else if (!strcmp(cmd,"execq"))
		{
			rc = GCALab_ProcessCommandQueue(cur_ws);
			GCALab_HandleErr(rc);
		}
		else if (!strcmp(cmd,"stopq"))
		{
			rc = GCALab_PauseCommandQueue(cur_ws);
			GCALab_HandleErr(rc);
		}
		else if(!strcmp(cmd,"printq"))
		{
			rc = GCALab_PrintCommandQueue(cur_ws);
			GCALab_HandleErr(rc);
		}
		else if(!strcmp(cmd,"quit"))
		{
			GCALab_ShutDown(GCALAB_SUCCESS);
		}
		else
		{
			fprintf(stderr,"Unknown Command [%s], type help for guidence\n",cmd);
		}

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

	opts[0] = GCALab_ParseCommandLineArgs(argc,argv);
	if (!(opts[0]))
	{
		return GCALAB_CL_PARSE_ERROR;
	}
	GCALab_numWS = 0;
	GCALab_mode = (opts[0])->mode;
	return GCALAB_SUCCESS;
}

/* GCALab_InteractiveMode(): Runs GCALab interactively
 */
void GCALab_GraphicsMode(GCALab_CL_Options* opts)
{
}

/* GCALab_BatchMode(): Runs Commands is batch mode
 * TODO: support a Job queue
 */
char GCALab_BatchMode(GCALab_CL_Options* opts)
{
}


/* GCALab_GetCommandCode(): converts the command string into a
 *                          command id
 */
unsigned int GCALab_GetCommandCode(char* cmd)
{
	if (!strcmp(cmd,"load"))
	{
		return GCALAB_LOAD;
	}
	else if(!strcmp(cmd,"save"))
	{
		return GCALAB_SAVE;
	}
	else if(!strcmp(cmd,"simulate"))
	{
		return GCALAB_SIMULATE;
	}
	else if (!strcmp(cmd,"gca"))
	{
		return GCALAB_GCA;
	}
	else if (!strcmp(cmd,"entropy"))
	{
		return GCALAB_ENTROPY;
	}
	else if (!strcmp(cmd,"param"))
	{
		return GCALAB_PARAM;
	}
	else if (!strcmp(cmd,"reverse"))
	{
		return GCALAB_REVERSE;
	}
	else
	{
		return GCALAB_NOP;
	}

}

/* GCALab_Worker(): Processing thread attached to a workspace
 */
void * GCALab_Worker(void *params)
{
	unsigned char ws_id;
	unsigned char exit,error;
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
				break;
			case GCALAB_WS_STATE_PROCESSING:
				/*the next command*/
				GCALab_DoNextCommand(ws_id);
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
	GCALabOutput *res;

	res = NULL;
	
	GCALab_LockWS(ws_id);
	index = WS(ws_id)->qhead;
	cmd_id = WS(ws_id)->commandqueue[index];
	trgt_id = WS(ws_id)->commandtarget[index];
	nparams = WS(ws_id)->numparams[index];
	params = strvncpy(WS(ws_id)->commandparams[index],nparams,GCALAB_MAX_STRLEN);
	GCALab_UnLockWS(ws_id);

	/*do processing*/
	switch(cmd)
	{
		case GCALAB_NOP:
			return GCALAB_SUCCESS;
			break;
		case GCALAB_LOAD:
		{
			char *filename;
			filename = NULL;
			for (i=0;i<nparams;i++)
			{
				if(!strcmp(params[i],"-f"))
				{
					filename = params[++i];
				}
			}
		}
			break;
		case GCALAB_SAVE:
		{
			char *filename;
			unsigned char gca_res_flag;
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
				else if (!strcmp(para,s[i],"-g"))
				{
					gca_res_flag = 1;
				}
			}
		}
			break;
		case GCALAB_SIMULATE:
		{
			unsigned int Tfinal;
			unsigned char reInit,ic_type;
			char *ic_filename;
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
		}
			break;
		case GCALAB_GCA:
		{
			char *meshfile;
			unsigned int NCells,genus,windowsize;
			unsigned char rule_type,rule_code,s,k,eca;
		}
			break;
		case GCALAB_ENTROPY:
		{
			unsigned int numSamples;
			unsigned int type;
		}
			break;
		case GCALAB_PARAM:
		{
			unsigned int type;
		}
			break;
		case GCALAB_REVERSE:
		{
			
		}
			break;
	}

	/*clean up copied params*/
	for (i=0;i<nparams;i++)
	{
		free(params[i]);
	}
	free(params);
	
	GCALab_LockWS(ws_id);
	WS(ws_id)->qhead = (WS(ws_id)->qhead+1)%GCALAB_COMMAND_BUFFER_SIZE;
	WS(ws_id)->numcommands--;
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
	for (i=0;i<GCALab_numWS;i++)
	{
		GCALab_SetState(i,GCALAB_WS_STATE_EXITING);
	}
	pthread_exit(0);
}

/* GCALab_PrintHelp(): Prints the help menu
 */
void GCALab_PrintHelp(void)
{
	fprintf(stdout,"\nCommand List:\n");
	fprintf(stdout,"-------------\n");
	fprintf(stdout,"Command\t\targs\tDescription\n");
	fprintf(stdout,"nwork\t\tn\tCreate a new workspace with enough slots to host n GCAs.\n");
	fprintf(stdout,"printwork\tnone\tPrint current workspace.\n");
	fprintf(stdout,"listallwork\tnone\tPrint summary of all workspaces.\n");
	fprintf(stdout,"help\t\tnone\tPrints this help menu.\n");
	fprintf(stdout,"chwork\t\tid\tChanges current workspace to be workspace id.\n");
	fprintf(stdout,"q\t\tcmd\tEnqueues the GCA operation cmd in the current command queue.\n");
	fprintf(stdout,"cancel\t\tcmd_id\tSets GCA operation cmd_id to be ignored.\n");
	fprintf(stdout,"stopq\t\tnone\tThe current queue will pause after completion of the current operation.\n");
	fprintf(stdout,"execq\t\tnone\tThe current queue will restart processing.\n");
	fprintf(stdout,"printq\t\tnone\tPrints current queue.\n");
	fprintf(stdout,"quit\t\tnone\tExits GCALab.\n");

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
	switch(GCALab_Global[ws_id]->state)
	{
		case GCALAB_WS_STATE_IDLE:
			fprintf(stdout,"Idle\n");
			break;
		case GCALAB_WS_STATE_PROCESSING:
			fprintf(stdout,"Running\n");
			break;
		case GCALAB_WS_STATE_PAUSED:
			fprintf(stdout,"Paused\n");
			break;
		case GCALAB_WS_STATE_EXITING:
			fprintf(stdout,"Quitting\n");
			break;
		case GCALAB_WS_STATE_ERROR:
			fprintf(stdout,"Error\n");
			break;
	}
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
		switch(GCALab_Global[i]->state)
		{
			case GCALAB_WS_STATE_IDLE:
				fprintf(stdout,"I\n");
				break;
			case GCALAB_WS_STATE_PROCESSING:
				fprintf(stdout,"R\n");
				break;
			case GCALAB_WS_STATE_PAUSED:
				fprintf(stdout,"P\n");
				break;
			case GCALAB_WS_STATE_EXITING:
				fprintf(stdout,"Q\n");
				break;
			case GCALAB_WS_STATE_ERROR:
				fprintf(stdout,"E\n");
				break;
		}
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
