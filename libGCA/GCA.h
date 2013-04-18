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
/* File: GCA.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 26/03/2012
 * Last Modified: 18/01/2013
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
#define LIFE_RULE_TYPE 3

#ifndef DEFAULT_RULE_TYPE
	#define DEFAULT_RULE_TYPE CODE_RULE_TYPE
#endif 

#define PACKED_STORAGE_TYPE 0
#define UNPACKED_STORAGE_TYPE 1

#ifndef DEFAULT_STORAGE_TYPE
	#define DEFAULT_STORAGE_TYPE PACKED_STORAGE_TYPE
#endif

#define S_MAX 256

#define DEFAULT_NEIGHBOURHOOD_SIZE 3
#define DEFAULT_NUM_CELLS 32

#define VON_NEUMANN_NEIGHBOURHOOD_TYPE 0
#define MOORE_NEIGHBOURHOOD_TYPE 1

#ifndef DEFAULT_NEIGHBOURHOOD_TYPE
	#define DEFAULT_NEIGHBOURHOOD_TYPE VON_NEUMANN_NEIGHBOURHOOD_TYPE
#endif

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
	unsigned int LUT_size;
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

unsigned int *GenerateTopology(unsigned int *N,unsigned char *k, unsigned char nh_type, mesh *m);
CellularAutomatonParameters *CreateCAParams(unsigned char nh_type,mesh *m,state s,unsigned char rule_type, unsigned int rule, unsigned int ws);
GraphCellularAutomaton *CreateECA(unsigned int N,unsigned int k,unsigned int rule,unsigned int ws);
GraphCellularAutomaton *CreateGCA(CellularAutomatonParameters *params);
GraphCellularAutomaton *CopyGCA(GraphCellularAutomaton *GCA);

/* cell and config get/sets functions*/
void SetCAIC(GraphCellularAutomaton *GCA,chunk *ic,unsigned char type);
void ResetCA(GraphCellularAutomaton *GCA);
state GetCellStatePacked(GraphCellularAutomaton *GCA, unsigned int i,unsigned int t);
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
unsigned char NhElim(GraphCellularAutomaton *GCA,unsigned char *flags,state *theta_i,state *theta_j,unsigned int startcell);
unsigned char *GetFlags(GraphCellularAutomaton *GCA);
unsigned char IsGOE(GraphCellularAutomaton *GCA);
unsigned char isValid(GraphCellularAutomaton *GCA,chunk *config,unsigned char *flags);

/*Analysis functions*/
float ShannonEntropy(GraphCellularAutomaton *GCA, unsigned int T,float *pm,float* logs_pm,float *S_im,unsigned char *TFm,unsigned int *cm);
float WordEntropy(GraphCellularAutomaton *GCA,unsigned int T, float *pm,float* logs_pm,float *W_im,unsigned char *TFm, unsigned int *cm,unsigned int *wlm);
unsigned int *SumCAImages(GraphCellularAutomaton *GCA,unsigned int *counts,chunk *preImages,unsigned int n);
float* ComputeExactProbs(GraphCellularAutomaton *GCA);
float* InputEntropy(GraphCellularAutomaton *GCA,unsigned int T,float* mu, float* sigma,unsigned int *Qm, float *logQm,float *IEm);
float lambda_param(GraphCellularAutomaton *GCA);
float Z_param(GraphCellularAutomaton *GCA);
float G_density(GraphCellularAutomaton *GCA,chunk* ics, unsigned int n);
float AttLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned t);
float TransLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned t);
float* PopDensity(GraphCellularAutomaton *GCA,chunk* ics,unsigned int T, float *dense);

#endif
