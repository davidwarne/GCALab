/* File: GCA.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 26/03/2012
 * Last Modified: 14/05/2012
 *
 * Descritpion: Definition of Graph Cellular Automata Libarary
 *
 * =============================================================================
 */

#ifndef __GCA_H
#define __GCA_H

/*standard headers*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*custom headers*/
#include "mesh.h"

#define THRESH_RULE_TYPE 0
#define COUNT_RULE_TYPE 1
#define CODE_RULE_TYPE 2 

#ifndef DEFAULT_RULE_TYPE
	#define DEFAULT_RULE_TYPE CODE_RULE_TYPE
#endif 

#define PACKED_STORAGE_TYPE 0
#define UNPACKED_STORAGE_TYPE 1

#ifndef DEFAULT_STORAGE_TYPE
	#define DEFAULT_STORAGE_TYPE PACKED_STORAGE_TYPE
#endif

#define S_MAX 256

#define DEFAULT_NEIGHBOURHOOD 3
#define DEFAULT_NUM_CELLS 32

#define POINT_IC_TYPE 0
#define NOISE_IC_TYPE 1
#define STRIPE_IC_TYPE 2
#define CHECKER_IC_TYPE 3
#define EXPLICIT_IC_TYPE 4

#define MAX_PRE_IMAGE_RETURN 1000

#ifndef DEFAULT_WINDOW_SIZE
	#define DEFAULT_WINDOW_SIZE 1200
#endif

#ifndef DEFAULT_IC_TYPE
	#define DEFAULT_IC_TYPE POINT_IC_TYPE
#endif

#ifndef CHUNK_SIZE_BITS
	#define CHUNK_SIZE_BITS 32
	typedef unsigned int chunk;
	//typedef unsigned long long chunk;
#else
	#if CHUNK_SIZE_BITS == 64
		typedef unsigned long long chunk;
	#elif CHUNK_SIZE_BITS == 32
		typedef unsigned int chunk;
	#elif CHUNK_SIZE_BITS == 16
		typedef unsigned short chunk;
	#elif CHUNK_SIZE_BITS == 8
		typedef unsigned char chunk;
	#endif
#endif

#define FACE_CELL_TYPE 0
#define VERT_CELL_TYPE 1
#ifndef DEFAULT_CELL_TYPE
	#define DEFAULT_CELL_TYPE VERT_CELL_TYPE
#endif

typedef unsigned char state;

typedef struct GraphCellularAutomaton_struct GraphCellularAutomaton;
typedef struct CellularAutomatonParameters_struct CellularAutomatonParameters; 

struct CellularAutomatonParameters_struct
{
	unsigned int N; /*number of cells*/
	unsigned int WSIZE; /*number of timesteps to store*/
	unsigned int rule; /* defines CA evolution rule*/
	state s;/*number of states/symbols*/
	unsigned char rule_type; /* defines type of rule*/
	unsigned char k;/*size of neigbourhood*/
	unsigned int *graph; /*topology of CA*/
};

/*Defintion of our graph CA*/
struct GraphCellularAutomaton_struct
{
	unsigned char log2s; /*number of bits per symbol*/
	unsigned char LUT_size;
	unsigned int size; /*number of chunks used to store the configuration*/
	unsigned int t; /*timestep*/
	state *ruleLUT; /*rule look-up table*/
	chunk *ic; /*the last initial condition set*/
	chunk *config; /*current configuration*/	
	chunk **st_pattern;/*spatio-temporal pattern*/
	CellularAutomatonParameters *params;
};

/*function prototypes*/

/*CA creation functions*/
unsigned int *GenerateTopology(unsigned int *N,unsigned char *k, mesh *m);
GraphCellularAutomaton *CreateECA(unsigned int N,unsigned int k,unsigned int rule);
GraphCellularAutomaton *CreateGCA(CellularAutomatonParameters *params);
GraphCellularAutomaton *CopyGCA(GraphCellularAutomaton *GCA);

/* cell and config get/sets functions*/
void SetCAIC(GraphCellularAutomaton *GCA,chunk *ic,unsigned char type);
void ResetCA(GraphCellularAutomaton *GCA);
state GetCellState(GraphCellularAutomaton *GCA, unsigned int i,unsigned int t);
void SetCellStatePacked(GraphCellularAutomaton *GCA, unsigned int i,state s);
state GetCellStatePacked_external(GraphCellularAutomaton *GCA,chunk* config, unsigned int i);
void SetCellStatePacked_external(GraphCellularAutomaton *GCA,chunk* config, unsigned int i,state s);
unsigned int* GetNeighbourhood(GraphCellularAutomaton * GCA,unsigned int i);
unsigned int GetNeighbourhood_config(GraphCellularAutomaton * GCA,unsigned int i,unsigned int t);
unsigned int GetNeighbourhood_config_external(GraphCellularAutomaton * GCA,chunk* config,unsigned int i);

/*Simulation functions*/
void CASimTSteps(GraphCellularAutomaton *GCA,unsigned int t);
unsigned int CANextStep(GraphCellularAutomaton *GCA);
chunk* CASimToAttCyc(GraphCellularAutomaton *GCA,unsigned int t);
unsigned char IsAttCyc(GraphCellularAutomaton *GCA);
chunk *CAGetPreImages(GraphCellularAutomaton *GCA,unsigned int* n,unsigned char* flags);
unsigned char *GetFlags(GraphCellularAutomaton *GCA);
unsigned char IsGOE(GraphCellularAutomaton *GCA);
unsigned char isValid(GraphCellularAutomaton *GCA,chunk *config,unsigned char *flags);

/*Analysis functions*/
float ShannonEntropy(GraphCellularAutomaton *GCA,unsigned int T, float *pr);
float WordEntropy(GraphCellularAutomaton *GCA,unsigned int T, float *pr);
unsigned int *SumCAImages(GraphCellularAutomaton *GCA,unsigned int *counts,chunk *preImages,unsigned int n);
float* ComputeExactProbs(GraphCellularAutomaton *GCA);
float* InputEntropy(GraphCellularAutomaton *GCA,unsigned int T,float *mu, float *sigma); 
float lambda_param(GraphCellularAutomaton *GCA);
float Z_param(GraphCellularAutomaton *GCA);
unsigned int G_density(GraphCellularAutomaton *GCA,chunk* ics, unsigned int n);

/*TODO: move this Thread stuff into another file, too lazy for now*/
#ifndef NO_THREADS
	#include <pthread.h>
	#ifndef NUM_THREADS
		#define NUM_THREADS 4
	#endif
	/*a struct to simply make passing data to threads easier*/
	struct ThreadData_struct
	{
		unsigned int ID;
		void * param1;
		void * param2;
		void * param3;
	};
	typedef struct ThreadData_struct ThreadData; 
#endif

#endif