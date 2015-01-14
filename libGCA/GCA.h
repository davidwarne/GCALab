/*
 * GCALab: An analysis tool for Graph Cellular Automata
 * Copyright (C) 2013  David J. Warne
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
/** 
 * @file GCA.h
 *
 * 
 *
 * @brief Definitions of Graph Cellular Automata Library.
 * @details The library contains functions for the creation and analysis of Graph Cellular 
 * Automata (i.e., Cellular Automata with a topology defined by a connected graph).
 * Analysis functions include entropy calculations, parameters such as Langton's lambda
 * and configuration transition graph metrics like the G-density.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Electrical Engineering and Computer Science
 * @author Faculty of Science and Engineering
 * @author Queensland University of Technology
 *
 * @version 0.19
 * @date 26/03/2012 - 18/01/2013
 * @copyright GNU Public License.
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

/** @brief Code to flag a Threshold type CA rule.*/
#define THRESH_RULE_TYPE 0
/** @brief Code to flag a totalistic type CA rule.*/
#define COUNT_RULE_TYPE 1
/** @brief Code to flag a Wolfram Code type CA rule.*/
#define CODE_RULE_TYPE 2 
/** @brief Code to flag a Game-of-life type CA rule (outer-totalisitc).*/
#define LIFE_RULE_TYPE 3

#ifndef DEFAULT_RULE_TYPE
/** @brief The rule type used by a CA if none is specified. */
	#define DEFAULT_RULE_TYPE CODE_RULE_TYPE
#endif 

/** @brief Flag to set CA memory usage to packed mode. */
#define PACKED_STORAGE_TYPE 0
/** @brief Flag to set CA memory usage to unpacked mode. */
#define UNPACKED_STORAGE_TYPE 1

#ifndef DEFAULT_STORAGE_TYPE
/** @brief The storage type used if none is specified. */
	#define DEFAULT_STORAGE_TYPE PACKED_STORAGE_TYPE
#endif

/** @brief The maximum number of states allowed for a CA.*/
#define S_MAX 256

/** @brief The CA neighbourhood size $k$ used when none is specified.*/
#define DEFAULT_NEIGHBOURHOOD_SIZE 3
/** @brief The CA cell count $n$ when none is specified.*/
#define DEFAULT_NUM_CELLS 32

/** @brief Code to flag a Von-Neumann CA neighbourhood type.*/
#define VON_NEUMANN_NEIGHBOURHOOD_TYPE 0
/** @brief Code to flag a Moore CA neighbourhood type.*/
#define MOORE_NEIGHBOURHOOD_TYPE 1

#ifndef DEFAULT_NEIGHBOURHOOD_TYPE
/** @brief The CA neighbourhood type used when none is specified.*/
	#define DEFAULT_NEIGHBOURHOOD_TYPE VON_NEUMANN_NEIGHBOURHOOD_TYPE
#endif

/** @brief Code to flag a single point initial condition.*/
#define POINT_IC_TYPE 0
/** @brief Code to flag a random noise initial condition.*/
#define NOISE_IC_TYPE 1
/** @brief Code to flag a striped initial condition.*/
#define STRIPE_IC_TYPE 2
/** @brief Code to flag a checker board initial condition.*/
#define CHECKER_IC_TYPE 3
/** @brief Code to flag a user defined initial condition.*/
#define EXPLICIT_IC_TYPE 4

#ifndef DEFAULT_IC_TYPE
/** @brief The initial condition type used when none is specified.*/
	#define DEFAULT_IC_TYPE POINT_IC_TYPE
#endif

/** @brief Limit on the number of pre-images returned by the EDEN-DET() algorithm.*/
#define MAX_PRE_IMAGE_RETURN 1000

#ifndef DEFAULT_WINDOW_SIZE
/** @brief The number of stored time steps if none is specified.*/
	#define DEFAULT_WINDOW_SIZE 1200
#endif


#ifndef CHUNK_SIZE_BITS
    /** @brief The number of bits in a memory chunk*/
	#define CHUNK_SIZE_BITS 32
    /** @brief A memory chunk used for CA state storage*/
	typedef unsigned int chunk;
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

/** @brief Code to flag that cells are derived from mesh faces.*/
#define FACE_CELL_TYPE 0
/** @brief Code to flag that cells are derived from mesh vertices.*/
#define VERT_CELL_TYPE 1
#ifndef DEFAULT_CELL_TYPE
/** @brief The cell derivation method used if none is specified.*/
	#define DEFAULT_CELL_TYPE VERT_CELL_TYPE
#endif

/** @brief A CA cell state.*/
typedef unsigned char state;

/** @brief A Graph Cellular Automaton.*/
typedef struct GraphCellularAutomaton_struct GraphCellularAutomaton;
/** @brief Graph Cellular Automaton Parameters.*/
typedef struct CellularAutomatonParameters_struct CellularAutomatonParameters; 

/** @brief A Graph Cellular Automaton parameter structure.*/
struct CellularAutomatonParameters_struct
{
	/** @brief Number of cells*/
	unsigned int N;
	/** @brief Number of timesteps to store*/
	unsigned int WSIZE; 
	/** @brief  Defines CA evolution rule*/
	unsigned int rule; 
	/** @brief Number of states/symbols*/
	state s;
	/** @brief  Defines type of rule*/
	unsigned char rule_type; 
	/** @brief Size of neigbourhood*/
	unsigned char k;
	/** @brief Topology of CA*/
	unsigned int *graph; 
};

/** @brief A Graph Cellular Automaton structure.*/
struct GraphCellularAutomaton_struct
{
	/** @brief Number of bits per symbol.*/
	unsigned char log2s; 
	/** @brief Size of state transition rule lookup table.*/
	unsigned int LUT_size;
	/** @brief Number of chunks used to store the configuration.*/
	unsigned int size; 
	/** @brief Current timestep.*/
	unsigned int t; 
	/** @brief Rule look-up table.*/
	state *ruleLUT; 
	/** @brief The last initial condition set.*/
	chunk *ic; 
	/** @brief The current configuration.*/	
	chunk *config; 
	/** @brief Spatio-temporal pattern.*/
	chunk **st_pattern;
	/** @brief Cellular Automaton Parameters.*/
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
void RotateNeighbourhood(GraphCellularAutomaton * GCA, unsigned int i, unsigned int r);
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
float AttLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned int t);
float TransLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned int t);
float* PopDensity(GraphCellularAutomaton *GCA,chunk* ics,unsigned int T, float *dense);

#endif
