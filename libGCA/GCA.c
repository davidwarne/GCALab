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
/* File: GCA.c
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 26/03/2012
 * Last Modified: 02/09/2012
 *
 * Version History:
 *       v 0.01 (26/03/2012) - i. Initial Version... feel like some hackin'
 *
 *       v 0.02 (26/03/2012) - i. Implemented GetCellStatePacked(), 
 *                                and SetCellStatePacked(). 
 *                             ii. Partially implemented GenerateTopology(), 
 *                                and CreateGCA()
 *
 *       v 0.03 (28/03/2012) - i. Implemented CreateECA(), CASimTSteps() 
 *                                 and CANextStep()
 *                             ii. did some testing of these functions.
 *
 *       v 0.04 (13/04/2012) - i. Implemented ShannonEntropy(),WordEntropy(),
 *                                and SumCAImages(). 
 *                             ii. added IC types to SetCAIC() 
 *                             iii. Completed TODO's 1,2, and 3. 
 *                             iv. Added 4,5,6, and 7 to TODO list
 *
 *       v 0.05 (16/04/2012) - i. Changed CANextStep() to work only with LUT, and
 *                                fixed CreateCA() to create LUT's for all rule defs.
 *                             ii. Implemented ComputeExactProbs() and 
 *                                 ComputeExactProbs_worker() using pthreads
 *
 *       v 0.06 (25/04/2012) - i. Implemented CAGetPreImages(), a really neat 
 *                                algorithm for finding the all the Pre-Images
 *                                directly. This is not based on Wuensche's method.
 *
 *       v 0.07 (14/05/2012) - i. Added more get/set configuration functions.
 *                             ii. broke up Getflags() and CAGetPreImages() as sometimes
 *                                 you don't actually want the Pre-images and are just
 *                                 testing for GOE.
 *                             iii. Added 4 to the Known Issue list
 *
 *       v 0.08 (15/05/2012) - i. Fixed problem with CAGetPreImages()
 *
 *       v 0.09 (19/06/2012) - i. Yay! finally get to do more coding, stupid project proposals...
 *                             ii. HPC upgrade is due to complete this week, so if I get my entropy 
 *                                 code ready I get free reign of all 1200 cores... sweet
 *                             iii. fixed SumImages(), it was useless for entropy evaluation, 
 *                                  this version will be much better.
 *                             iv. added IsAttCyc() and IsGOE() functions, just to avoid unnecessar 
 *                                 overheads of using CAGetPreImages().
 *
 *       v 0.10 (1/09/2012) - i. Found unfortunate issue with my IsGOE method, added issue 5.
 *                            ii. Implemented lambda_param() and InputEntropy().
 *
 *       v 0.11 (2/09/2012) - i. Implemented G_density() and Z_param().
 *       
 *       v 0.12 (14/10/2012) - i. separated RevAlgIter() from the GetFlags() function
 *                             ii. worked on a solution to Known Issue 5, the idea
 *                                 is based on the fact that nhoods that cause in-
 *                                 consistencies are eliminated in one iteration of
 *                                 RevAlgIter() if they are assumed. This helps the 
 *                                 entire algorithm to converge such that 
 *                                 rowOfZeros <=> GOE. Considering this issue resolved.
 *                             iii. I have been considering names for my rev algorithm
 *                                  , maybe neighbourhood elimination?
 *       v 0.13 (19/10/2012) - i. fixed bug in Word Entropy causing seg fault.
 *       
 *       v 0.14 (20/10/2012) - i. fixed bug in Input Entropy.
 *                             ii. fixed memory errors in CreateGCA.
 *
 * Description: Implementation of Graph Cellular Automata Libarary
 *
 * TODO List:
 *    1. Implement the m != NULL case for GenerateTopology(). - done (v 0.04)
 *    2. Add other rule type generators in CreateGCA(). - done (v 0.04)
 *    3. Update CANextStep() for other rule types. - done (v 0.04)
 *    4. SetCAIC is a bit ECA centric, should extend to usefule ICs for graphs
 *       in general.
 *    5. ShannonEntropy() needs to be rigorously tested - done (v 0.05).
 *    6. WordEntropy() needs to be rigorously tested.
 *    7. The CANextStep() function should be benchmarked for different rule types,
 *       if no advantage in reduced formats then an explicit neighbourhood 
 *       lookup may as well be done to simplify this function - done (v 0.05)
 *    8. Implement my reverse algorithm, should be faster than Wuensche's and you
 *       get all the pre-images for free - done (v 0.06).
 *
 * Known Issues:
 *     1. Currently code geared toward space optimisation, lots of bit twiddling
 *        this may be a performance issue later. I have Pre-processor defintions 
 *        for packed and unpacked types but only packed are implemented.
 *     2. Threshold and Count type rules are limited to a single state of 
 *        influence, only a problem for non-binary CA.
 *     3. Calculate Exact Probabilities is wrong... a whole new methodology 
 *        is required - fixed (v 0.09)
 *     4. There is a bug in CAGetPreImages(), it is when calculating permutations
 *        of possible preimages. - fixed (v 0.08)
 *     5. I've managed to prove by counter-example (sadly) that my GOE test condition (
 *        column of zeros) is a sufficient but not necessary condition for the existence.
 *        That is colOfzeros => GOE not colOfzeros <=> GOE... working on a further condition
 *        such that colOfzeros ^ q <=> GOE. - fixed (v 0.12)
 *     6. There seems to be a bug in the InputEntropy Calculation, need to look into this.
 *     7. Word entropy is causing seg faults. - fixed (v 0.13)
 * =============================================================================
 */

#include "GCA.h"

/* GenerateTopology(): Creates a topology array from a mesh. If meash is NULL
 *                     then a regular 1-dimensional genus-1 topology of N cells
 *                     is created.
 * Parameters:
 *     m - mesh structure to generate the topology from, can be NULL
 *     N - pointer to number of cells
 *     k - pointer to neighbourhood size
 * Returns:
 *     an array of unsigned int of length m->vList->numVerts)*k or N*k
 */
unsigned int * GenerateTopology(unsigned int *N, unsigned char *k,mesh *m)
{
	unsigned int *graph;
	unsigned int i,j,ii,jj,r,*e;
	unsigned int N_tmp,k_tmp;
	int *f1,*f2;
	unsigned char tf;
	if (m == NULL)
	{
		N_tmp = (N == NULL) ? DEFAULT_NUM_CELLS : *N;
		k_tmp = (k == NULL) ? DEFAULT_NEIGHBOURHOOD : *k;
		
		graph = (unsigned int *)malloc(N_tmp*(k_tmp-1)*sizeof(unsigned int));
		if (!graph)
		{
			return NULL;
		}
		memset((void*)graph,0,N_tmp*(k_tmp-1)*sizeof(unsigned int));
		
		/*neighbourhood radius */
		r = (k_tmp-1)/2;
		
		/*periodic boundaries*/
		for (i=0;i<N_tmp;i++)
		{
			for (j=0;j<r;j++)
			{
				graph[i*(k_tmp-1) + j] = (i + j - r + N_tmp)%N_tmp;
			}
			for (j=r;j<2*r;j++)
			{
				graph[i*(k_tmp-1) + j] = (i + j-r+1)%N_tmp;
			}
		}
	}
	else
	{
		N_tmp = m->fList->numFaces;
		k_tmp = m->fList->maxVerts+1;
		graph = (unsigned int *)malloc(N_tmp*(k_tmp-1)*sizeof(unsigned int));
		if (!graph)
		{
			return NULL;
		}
		memset((void*)graph,0xFF,N_tmp*(k_tmp-1)*sizeof(unsigned int));
		
		/*neighbours are the adjacent faces*/
		for (i=0;i<N_tmp;i++)
		{
			f1 = GetFace_ptr(m->fList,i);
			for (j=0;j<k_tmp-2;j++)
			{
				for (ii=0;ii<N_tmp;ii++)
				{
					f2 = GetFace_ptr(m->fList,ii);
					tf = 0;
					for (jj=0;jj<k_tmp-2;jj++)
					{
						if (((f1[j] == f2[jj]) && (f1[j+1] == f2[jj+1])) 
							|| ((f1[j] == f2[jj+1]) && (f1[j+1] == f2[jj])))
						{
							tf |= 1;
						}
					}
					if (((f1[j] == f2[k_tmp-2]) && (f1[j+1] == f2[0])) 
						|| ((f1[j] == f2[0]) && (f1[j+1] == f2[k_tmp-2])))
					{
						tf |= 1;
					}

					if (tf && (ii != i))
						break;
				}	
				graph[i*(k_tmp-1) + j] = ii;
			}
			
			for (ii=0;ii<N_tmp;ii++)
			{
				f2 = GetFace_ptr(m->fList,ii);
				tf = 0;
				for (jj=0;jj<k_tmp-2;jj++)
				{
					if (((f1[k_tmp-2] == f2[jj]) && (f1[0] == f2[jj+1])) 
						|| ((f1[k_tmp-2] == f2[jj+1]) && (f1[0] == f2[jj])))
					{
						tf |= 1;
					}
				}
				if (((f1[k_tmp-2] == f2[k_tmp-2]) && (f1[0] == f2[0])) 
					|| ((f1[k_tmp-2] == f2[0]) && (f1[0] == f2[k_tmp-2])))
				{
					tf |= 1;
				}
				if (tf && (ii != i))
					break;
			}	
			graph[i*(k_tmp-1) + k_tmp-2] = ii;
		}
	}
	
	*N = N_tmp;
	*k = k_tmp;
	return graph;
}

/* CreateCAParams(): Creates consistent parameters for a Graph Cellular Automata on a given 
 *                   topology.
 * Parameters:
 *     m - the Geometry of the domain (which is used to derive topology)
 *     rule_type - the type of rule that the rule code represents
 *     rule - the rule code
 *     ws - the window-size (i.e., the number of timesteps to store)
 */
CellularAutomatonParameters *CreateCAParams(mesh *m,state s,unsigned char rule_type, unsigned char rule, unsigned int ws)
{
	CellularAutomatonParameters *params;

	params = (CellularAutomatonParameters *)malloc(sizeof(CellularAutomatonParameters));
	if (!params)
	{
		return NULL;
	}

	params->graph = GenerateTopology(&(params->N),&(params->k),m);
	if (!(params->graph))
	{
		return NULL;
	}

	params->WSIZE = (ws == 0) ? DEFAULT_WINDOW_SIZE : ws;
	params->rule_type = rule_type;
	params->rule = rule;
	params->s = s;
	return params;
}

/* CreateECA(): Creates a Graph representation of an Elementary Cellular 
 *              Automaton
 *
 * Parameters:
 *     N - number of cells
 *     k - neighbour size (i.e. 2*r +1 )
 *     rule - Wolfram's ECA rule notation
 *
 * Returns:
 *     A Graph Elementary Cellular Automaton ready for simlulation
 */
GraphCellularAutomaton *CreateECA(unsigned int N,unsigned int k,unsigned int rule)
{
	CellularAutomatonParameters *params;
	GraphCellularAutomaton *ECA;
	unsigned int N_tmp;
	unsigned char k_tmp;
	
	N_tmp = N;
	k_tmp = (unsigned char)k;
	params = (CellularAutomatonParameters *)malloc(sizeof(CellularAutomatonParameters));
	if (!params)
	{
		return NULL;
	}
	
	params->s = 2; /*binary {0,1}*/
	params->k = (unsigned char)k;
	params->N = N;
	params->WSIZE = DEFAULT_WINDOW_SIZE;
	params->rule = rule;
	params->rule_type = CODE_RULE_TYPE;
	params->graph = GenerateTopology(&N_tmp,&k_tmp,NULL);
	
	ECA = CreateGCA(params);
	
	return ECA;
}

/* CreateGCA(): Creates a Graph Cellular Automaton
 *
 * Parameters:
 *     params - parameters defining CA dynamics and topology
 *
 * Returns:
 *     A Graph Cellular Automaton ready for simlulation
 */
GraphCellularAutomaton *CreateGCA(CellularAutomatonParameters *params)
{
	GraphCellularAutomaton *GCA;
	unsigned int i,j;
	/*allocate memory*/
	if (!(GCA = (GraphCellularAutomaton *)malloc(sizeof(GraphCellularAutomaton))))
	{
		return NULL;
	}
	
	
	GCA->log2s = 0;
	GCA->params = params;
	/* calculate the number of bits per symbol (needed alot later)*/
	{ 
		register state s; s = params->s;
		while (s >>= 1) GCA->log2s++;
	}
	
	GCA->size = ceil((float)(((params->N)*(GCA->log2s))) / (float)CHUNK_SIZE_BITS);
	if (!(GCA->st_pattern = (chunk **)malloc((params->WSIZE)*sizeof(chunk *))))
	{
		return NULL;
	}
	
	for (i=0;i<(params->WSIZE);i++)
	{
		GCA->st_pattern[i] = (chunk *)malloc((GCA->size)*sizeof(chunk));
		if (!(GCA->st_pattern[i]))
		{
			return NULL;
		}
		for (j=0;j<(GCA->size);j++)
		{
			GCA->st_pattern[i][j] = 0;
		}
	}
	
	GCA->config = GCA->st_pattern[0];
	GCA->t = 0;
	/*intial condition copy memory*/
	GCA->ic = (chunk *)malloc((GCA->size)*sizeof(chunk));
	if (!(GCA->ic))
	{
		return NULL;
	}
	/*create the rule table*/
	switch(params->rule_type)
	{
		case THRESH_RULE_TYPE: /*for now based on a single state*/
		{	
			state ref;
			unsigned char thresh;
			GCA->LUT_size = pow(params->s,params->k);
			GCA->ruleLUT = (state *)malloc((GCA->LUT_size)*sizeof(state));
			
			if (!(GCA->ruleLUT))
			{
				return NULL;
			}
			/* here the rule is composed of 5 states-- the state of interest, 
			 * the threshhold,follow by three states <,=,>
			 */ 
			ref = (state)(GCA->params->rule & ((0x1 << (GCA->log2s)) - 1));
			thresh = (unsigned char)((GCA->params->rule >> (GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1)); 
			for (i=0;i<GCA->LUT_size;i++)
			{
				register unsigned char sum;
				register state s_i;
				register unsigned char ii;
				sum = 0;
				for (j=0;j<GCA->params->k;j++)
				{
					s_i = (state)((i >> i*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
					sum += (s_i == ref);
				}
				
				if (sum < thresh)
				{
					ii = 2;
				}
				else if (sum == thresh)
				{
					ii = 3;
				}
				else
				{
					ii = 4;
				}
				GCA->ruleLUT[i] = (state)((GCA->params->rule >> ii*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
			}
		}
			break;
		case COUNT_RULE_TYPE:
		{
			state ref;
			GCA->LUT_size = pow(params->s,params->k);
			GCA->ruleLUT = (state *)malloc((GCA->LUT_size)*sizeof(state));
		
			if (!(GCA->ruleLUT))
			{
				return NULL;
			}
			
			/* here the rule is composed of k+1 states-- the state of interest, 
			 * follow by k state changes for each count.
			 */ 
			ref = (state)(GCA->params->rule & ((0x1 << (GCA->log2s)) - 1));
			
			for (i=0;i<GCA->LUT_size;i++)
			{
				register unsigned char sum;
				register state s_i;
				register unsigned char ii;
				sum = 0;
				for (j=0;j<GCA->params->k;j++)
				{
					s_i = (state)((i >> i*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
					sum += (s_i == ref);
				}
				
				for (j=0;j<GCA->params->k;j++)
				{
					if (sum == (state)j)
					{
						ii = (unsigned char)j;
						break;
					}
				}
				
				GCA->ruleLUT[i] = (state)((GCA->params->rule >> ii*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
			}
		}
			break;
		case CODE_RULE_TYPE:/*rule is the Wolfram rule code*/
			GCA->LUT_size = pow(params->s,params->k);
			GCA->ruleLUT = (state *)malloc((GCA->LUT_size)*sizeof(state));

			if (!(GCA->ruleLUT))
			{
				return NULL;
			}
			
			for (i=0;i<GCA->LUT_size;i++)
			{
				GCA->ruleLUT[i] = params->rule >> i*(GCA->log2s);		
				GCA->ruleLUT[i] &= ((0x1 << (GCA->log2s)) - 1); 
			}
			break;
	}
	
	return GCA;
}

/* CopyGCA(): Creates copy of the given Graph Cellular Automaton
 *
 * Parameters:
 *     GCA - The Graph Cellular Automaton to copy
 *
 * Returns:
 *     A Graph Cellular Automaton identical to the given one
 */
GraphCellularAutomaton *CopyGCA(GraphCellularAutomaton *GCA)
{
	GraphCellularAutomaton *GCA_cp;
	int i,j;
	/*easy case... :) */
	if (GCA == NULL)
	{
		return NULL;
	}
	
	GCA_cp = (GraphCellularAutomaton *)malloc(sizeof(GraphCellularAutomaton));
	if (!GCA_cp)
	{
		return NULL;
	}
	
	GCA_cp->params = (CellularAutomatonParameters *)malloc(sizeof(CellularAutomatonParameters));
	if (!(GCA_cp->params))
	{
		return NULL;
	}
	
	/*copy parameters*/
	GCA_cp->params->N = GCA->params->N;
	GCA_cp->params->WSIZE = GCA->params->WSIZE;
	GCA_cp->params->rule = GCA->params->rule;
	GCA_cp->params->s = GCA->params->s;
	GCA_cp->params->rule_type = GCA->params->rule_type;
	GCA_cp->params->k = GCA->params->k;
	GCA_cp->params->graph = (unsigned int *)malloc((GCA->params->N)*(GCA->params->k-1)*sizeof(unsigned int));
	if (!(GCA_cp->params->graph))
	{
		return NULL;
	}
	for (i=0;i<(GCA->params->N)*(GCA->params->k-1);i++)
	{
		GCA_cp->params->graph[i] = GCA->params->graph[i];
	}
	
	/*copy the CA now*/
	GCA_cp->log2s = GCA->log2s;
	GCA_cp->LUT_size = GCA->LUT_size;
	GCA_cp->size = GCA->size;
	GCA_cp->t = GCA->t;
	GCA_cp->ruleLUT = (state*)malloc((GCA->LUT_size)*sizeof(state));
	if (!(GCA_cp->ruleLUT))
	{
		return NULL;
	}
	for (i=0;i<GCA->LUT_size;i++)
	{
		GCA_cp->ruleLUT[i] = GCA->ruleLUT[i];
	}
	GCA_cp->st_pattern = (chunk**)malloc((GCA->params->WSIZE)*sizeof(chunk*));
	if (!(GCA_cp->st_pattern))
	{
		return NULL;
	}
	for (i=0;i<GCA->params->WSIZE;i++)
	{
		GCA_cp->st_pattern[i] = (chunk*)malloc((GCA->size)*sizeof(chunk));
		if (!(GCA_cp->st_pattern[i]))
		{
			return NULL;
		}
	}
	for (i=0;i<GCA->params->WSIZE;i++)
	{
		for(j=0;j<GCA->size;j++)
		{
			GCA_cp->st_pattern[i][j] = GCA->st_pattern[i][j];
		}
	}
	GCA_cp->config = GCA_cp->st_pattern[0];

	GCA_cp->ic = (chunk*)malloc((GCA->size)*sizeof(chunk));
	if(!(GCA_cp->ic))
	{
		return NULL;
	}
	for (i=0;i<GCA->size;i++)
	{
		GCA_cp->ic[i] = GCA->ic[i];
	}

	return GCA_cp;
}

/* SetCAIC(): Set Cellular Automaton intitial configuration. if ic == NULL
 *            then initial conditions of the given type are generated
 *
 * Parameters:
 *     GCA - Graph Cellular Automaton
 *
 * TODO: improve the IC types right now very ECA centric
 *
 */
void SetCAIC(GraphCellularAutomaton *GCA,chunk *ic,unsigned char type)
{
	unsigned int i;
	
	if (ic != NULL)
	{
		for (i=0;i<GCA->size;i++)
		{
			GCA->config[i] = ic[i];
		}
	}
	else
	{
		switch(type)
		{
			default:
			case POINT_IC_TYPE:
				GCA->config[0] = 0x1;
				for (i=1;i<GCA->size;i++)
				{
					GCA->config[i] = 0x0;
				}
				break;
			case NOISE_IC_TYPE:
				for (i=0;i<GCA->size;i++)
				{
					GCA->config[i] = rand();
				}
				break;
			case STRIPE_IC_TYPE:
				for(i=0;i<GCA->size;i++)
				{
					GCA->config[i] = 0x1;
				}
				break;
			case CHECKER_IC_TYPE:
				GCA->config[0] = 0x0;
				for(i=1;i<GCA->size;i++)
				{
					GCA->config[i] = ~GCA->config[i-1];
				}
				break;
		}
	}
	/*store a copy of the initial condition for restart*/
	for (i=0;i<GCA->size;i++)
	{
		GCA->ic[i] = GCA->config[i];
	}
}

/* ResetCA(): Resets Time-evolution for the Cellular Automaton to t0
 *            with the last used initial conditions.
 * Parameters:
 *     GCA - The Graph Cellular Automaton to reset
 */
void ResetCA(GraphCellularAutomaton *GCA)
{
	/*reset initial conditions*/
	SetCAIC(GCA,GCA->ic,EXPLICIT_IC_TYPE);
	/*reset time to 0*/
	GCA->t = 0;
}

/* GetCellStatePacked(): gets the state of the ith cell in the current
 *                       configuration
 */ 
state GetCellStatePacked(GraphCellularAutomaton *GCA, unsigned int i,unsigned int t)
{
	register unsigned int log2s,p,r,q;
	log2s = GCA->log2s;
	p = CHUNK_SIZE_BITS/log2s;
	r = (i%p)*log2s;
	q = i/p;
	/* what the crap? gotta love bit twiddling*/
	return (GCA->st_pattern[t][q] >> r) & ((0x1 << log2s) - 1);
}

/* SetCellStatePacked(): sets the state of the ith cell in the current
 *                       configuration to the given symbol
 */ 
void SetCellStatePacked(GraphCellularAutomaton *GCA, unsigned int i,state s)
{
	register unsigned int log2s,p,r,q;
	register chunk mask;
	
	mask = 0;
	log2s = GCA->log2s;
	p = CHUNK_SIZE_BITS/log2s;
	r = (i%p)*log2s;
	q = i/p;
	
	/* In a year from now you'll be thinking "What the hell was I thinking?!"*/
	mask = s << r;
	GCA->config[q] &= ~((((0x1 << log2s) - 1) << r));
	GCA->config[q] |= mask;
}

/* GetCellStatePacked_external(): gets the state of the ith cell in the given
 *                       configuration
 */ 
state GetCellStatePacked_external(GraphCellularAutomaton *GCA,chunk* config, unsigned int i)
{
	register unsigned int log2s,p,r,q;
	log2s = GCA->log2s;
	p = CHUNK_SIZE_BITS/log2s;
	r = (i%p)*log2s;
	q = i/p;
	/* what the crap? gotta love bit twiddling*/
	return (config[q] >> r) & ((0x1 << log2s) - 1);
}

/* SetCellStatePacked_external(): sets the state of the ith cell in the current
 *                       configuration to the given symbol will not update the GCA 
 */
void SetCellStatePacked_external(GraphCellularAutomaton *GCA,chunk* config, unsigned int i,state s)
{
	register unsigned int log2s,p,r,q;
	register chunk mask;
	mask = 0;
	log2s = GCA->log2s;
	p = CHUNK_SIZE_BITS/log2s;
	r = (i%p)*log2s;
	q = i/p;
	
	/* what the crap?*/
	mask = s << r;
	config[q] &= ~((((0x1 << log2s) - 1) << r));
	config[q] |= mask;
}

/* GetNeighbourhood(): gets the list of cells that are in the neighbourhood of
 *                     cell i
 *
 * Note: 
 *   returns a reference to the GCA graph.
 */
unsigned int* GetNeighbourhood(GraphCellularAutomaton * GCA,unsigned int i)
{
	return GCA->params->graph + i*(GCA->params->k-1);
}

/* GetNeighbourhood_config_external(): gets configuration of the given neighbourhood
 *
 */
unsigned int GetNeighbourhood_config_external(GraphCellularAutomaton * GCA,chunk* config,unsigned int i)
{
	unsigned int *U_i;
	unsigned int nhood,j;
	U_i  = GCA->params->graph + i*(GCA->params->k-1);
	nhood = 0;
	nhood |= GetCellStatePacked_external(GCA,config,i) << GCA->log2s*((GCA->params->k-1)/2);
		
	for (j=0;j<(GCA->params->k-1)/2;j++)
	{
		nhood |= GetCellStatePacked_external(GCA,config,U_i[j]) << GCA->log2s*j;
		
	}
	for (j=(GCA->params->k-1)/2;j<(GCA->params->k-1);j++)
	{
		nhood |= GetCellStatePacked_external(GCA,config,U_i[j]) << GCA->log2s*(j+1);
	}
	return nhood;
}

/* GetNeighbourhood_config(): gets configuration of the given neighbourhood
 *
 */
unsigned int GetNeighbourhood_config(GraphCellularAutomaton * GCA,unsigned int i,unsigned int t)
{
	unsigned int * U_i;
	unsigned int nhood,j;
	U_i  = GCA->params->graph + i*(GCA->params->k-1);
	nhood = 0;
	nhood |= GetCellStatePacked(GCA,i,t) << GCA->log2s*((GCA->params->k-1)/2);
		
	for (j=0;j<(GCA->params->k-1)/2;j++)
	{
		nhood |= GetCellStatePacked(GCA,U_i[j],t) << GCA->log2s*j;
		
	}
	for (j=(GCA->params->k-1)/2;j<(GCA->params->k-1);j++)
	{
		nhood |= GetCellStatePacked(GCA,U_i[j],t) << GCA->log2s*(j+1);
	}
	return nhood;
}

/* CASimTSteps(): Evolves a Cellular Automaton up to timestep t
 *
 * Parameters:
 *      GCA - a pointer to the graph Cellular Automaton
 *      t - timestep to terminate at.
 */
void CASimTSteps(GraphCellularAutomaton *GCA,unsigned int t)
{
	if (t < GCA->t)
	{
		ResetCA(GCA);
	}
	while(CANextStep(GCA) != t);
}

/* CANextStep(): Evolves the Cellular Automaton 1 timestep
 *
 * Parameters:
 *      GCA - a pointer to the graph Cellular Automaton
 */
unsigned int CANextStep(GraphCellularAutomaton *GCA)
{
	unsigned int i,j,N,k,WSIZE;
	chunk *next_config;

	N = GCA->params->N;
	k = GCA->params->k;
	WSIZE = GCA->params->WSIZE;
	/*update the window*/
	next_config = GCA->st_pattern[WSIZE-1];
	
	//GCA->st_pattern[1] = GCA->st_pattern[0];
	for (i=WSIZE-1;i>0;i--)
	{
		GCA->st_pattern[i] = GCA->st_pattern[i-1];
	}
	
	GCA->st_pattern[0] = next_config;
	GCA->config = GCA->st_pattern[0];
	
	for (i=0;i<N;i++)
	{
		register unsigned int nhood;
		nhood = GetNeighbourhood_config(GCA,i,1);
		SetCellStatePacked(GCA,i,GCA->ruleLUT[nhood]);
	}
	GCA->t++;
	
	return GCA->t;
}

/* CASimToAttCyc(): Simulates the CA until an attractor cycle begins.
 * 
 * Paramters:
 *	   GCA - graph cellular automaton
 *     t - max time step, abort if a cycle is not reached before this step
 *
 * Returns:
 *    The first configuration to be repeated, hence the start of the cycle
 */
chunk* CASimToAttCyc(GraphCellularAutomaton *GCA,unsigned int t)
{
	chunk *cycle;
	unsigned int i;
	cycle = (chunk*)malloc((GCA->size)*sizeof(chunk));
	if (!cycle)
	{
		return NULL;
	}

	while (!IsAttCyc(GCA)) 
	{
		CANextStep(GCA);
	}

	for (i=0;i<GCA->size;i++)
	{
		cycle[i] = GCA->config[i];
	}
	return cycle;
}

/* IsAttCyc(): detects if the CA has entered an attractor cycle.
 *
 * Paramaters:
 *     GCA - yes, I'll just assume you know what this is...
 * Returns:
 *     1 if an attractor cycle has begun, 0 otherwise 
 */
unsigned char IsAttCyc(GraphCellularAutomaton *GCA)
{
	unsigned int WSIZE,nbytes,t;
	WSIZE = GCA->params->WSIZE;
	nbytes = (GCA->size)*sizeof(chunk);
	
	for (t=GCA->t;t>0;t--)
	{
		/*if we have seen this configuation before, then we have entered an attractor cycle */
		if (!memcmp((void*)(GCA->config),(void*)(GCA->st_pattern[t]),nbytes)){
			/*woah!? dejavu... was it the same cat?*/
			return 1;
		}
	}
	return 0;
}

/* CAGetPreImages(): Gets all the Pre-Images that can lead to this CAs current
 *                   configuration.
 * Paramters:
 *    GCA - graph cellular automaton
 *
 * Returns:
 *    An array of configurations that are the pre-images of GCA's current configuration
 *    NULL if configuration is a Garden-of-Eden configuration
 * TODO: A bit inefficient, and we miss some valid pre-images sometimes
 */
chunk *CAGetPreImages(GraphCellularAutomaton *GCA,unsigned int* n,unsigned char* flags)
{
	chunk *preImages;
	unsigned int upper_bound;
	unsigned int numPreImages;
	unsigned char p_states[S_MAX];
	unsigned char* counts;
	unsigned char* cur;
	unsigned int carry;
	int i,j,k;
	unsigned int sum;
	state st;
	/*first generate flages, if we need to*/
	if (!flags)
	{
		flags = GetFlags(GCA);
		if (flags == NULL)
		{
			return NULL;
		}
	}
	
	counts = (unsigned char*)malloc(GCA->params->N*sizeof(unsigned char));
	cur = (unsigned char*)malloc(GCA->params->N*sizeof(unsigned char));
	memset((void*)cur,0,GCA->params->N*sizeof(unsigned char));
	
	/*check if Garden-of-Eden configuration, get absolute upped bound on number of Pre-images*/
	upper_bound = 1;
	for (i=0;i<GCA->params->N;i++)
	{
		
		sum = 0;
		for (j=0;j<GCA->params->s;j++)
		{
			p_states[j] = 0;
		}
		
		for (j=0;j<GCA->LUT_size;j++)
		{
			sum += flags[i*GCA->LUT_size + j];
			p_states[(j >> GCA->log2s*((GCA->params->k-1)/2)) & (0x1 << (GCA->log2s)) - 1] |= flags[i*GCA->LUT_size + j];
		}
		
		if (sum != 0)
		{
			sum = 0;
			for (j=0;j<GCA->params->s;j++)
			{
				sum += p_states[j];
			}
			counts[i] = sum;
			upper_bound *= sum;
		}
		else
		{
			/*a cell has no valid pre-nieghbourhoods -> GOE*/
			*n = 0;
			return NULL;
		}
	}
	
	if (upper_bound > MAX_PRE_IMAGE_RETURN)
	{
		fprintf(stderr,"Warning! Possibly large number of pre-images (%d), may not be able to return them all...\n",upper_bound);
		upper_bound = MAX_PRE_IMAGE_RETURN;
	}
	
	/*allocate memory for pre-images*/
	preImages = (chunk *)malloc(GCA->params->N*upper_bound*sizeof(chunk));
	numPreImages = 0;
	carry = 1;
	
	for (i=0;i<upper_bound;i++)
	{
		/*get the next Permutation*/
		carry = 1;
		for (j=0;j<GCA->params->N;j++)
		{
			if (carry + cur[j] == counts[j])
			{
				cur[j] = 0;
				carry = 1;
			}
			else
			{
				cur[j] += carry;
				break;
			}
		}
		
		for (j=0;j<GCA->params->N;j++)
		{
			st = 0;
			sum = 0;
			for (k=0;k<GCA->LUT_size && sum <= cur[j];k++)
			{
				sum += flags[j*GCA->LUT_size+k];
				st = k >> GCA->log2s*((GCA->params->k - 1)/2) & ((0x1 << (GCA->log2s)) - 1);
			}
			SetCellStatePacked_external(GCA,preImages+numPreImages*(GCA->size),j,st);
		}
		
		/*check it is valid*/
		if (isValid(GCA,preImages+numPreImages*(GCA->size),flags))
		{
			numPreImages++;
		}
	}
	free(counts);
	free(cur);
	*n = numPreImages;
	return preImages;
}

unsigned char RevAlgIter(GraphCellularAutomaton *GCA,unsigned char *flags,state *theta_i,state *theta_j,unsigned int startcell)
{
	unsigned char r;
	state *W_i,*W_j;
	unsigned int i,j,ii,jj;
	unsigned int q,p,pp;
	unsigned int sum,prev_sum;
	unsigned int *U_i,*U_j;
	unsigned int state_mask,mask_i,mask_j,mask_ii,mask_jj;
	
	state_mask = (0x1 << (GCA->log2s)) - 1;
	r = GCA->log2s*(GCA->params->k-1)/2;
	/*for every connection*/
	for (i=startcell+1;i<GCA->params->N;i++)
	{

		U_i = GetNeighbourhood(GCA,i);
		for (j=0;j<(GCA->params->k-1);j++)
		{
			U_j = GetNeighbourhood(GCA,U_i[j]);
			/*get the neighbour number for each cell in the other neighbourhood*/
			ii = 0;
			while (U_j[ii] != i) ii++;
			jj = 0;
			while (U_i[jj] != U_i[j]) jj++;
			ii = (ii >= (GCA->params->k-1)/2) ? ii+1 : ii;
			jj = (jj >= (GCA->params->k-1)/2) ? jj+1 : jj;
			/*create masks for selecting respective Pre-neighbourhood states*/
			mask_i = (state_mask << r); 
			mask_jj = (state_mask << GCA->log2s*jj);
			mask_j = (state_mask << r); 
			mask_ii = (state_mask << GCA->log2s*ii);
				
			/*create list for sub-configs for match test*/
			memset((void*)theta_i,0,(GCA->params->s)*(GCA->params->s)*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/
				if( flags[i*(GCA->LUT_size) + q])
				{
					p = (((q & mask_i) >> r) << GCA->log2s);
					p |= ((q & mask_jj) >> GCA->log2s*jj); 
					/*only set if this still a possible pre-neighbourhood*/
					theta_i[p] = 1;
				}
			}
				
			memset((void*)theta_j,0,(GCA->params->s)*(GCA->params->s)*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/						
				if (flags[U_i[j]*(GCA->LUT_size) + q])
				{
					p = (((q & mask_ii) >> GCA->log2s*ii) << GCA->log2s);
					p |= ((q & mask_j) >> r); 
					/*only set if this still a possible pre-neighbourhood*/
					theta_j[p] = 1;
				}
			}
				
			/* if theta_i not in theta_j then matching Pre-neighbourhoods 
			 * cannot contribute Pre-Image around cell_i, likewise for 
			 * theta_j not in theta_i 
			 * TODO: could repace theta_i and theta_j with bit strings
			 */
			for (p=0;p<(GCA->params->s)*(GCA->params->s);p++)
			{
			 	if (theta_i[p] && !theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_i*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						if (flags[i*(GCA->LUT_size) + q])
						{
							pp = (((q & mask_i) >> r) << GCA->log2s);
							pp |= ((q & mask_jj) >> GCA->log2s*jj); 
							flags[i*(GCA->LUT_size) + q] = (p == pp) ? 0 : flags[i*(GCA->LUT_size) + q];
						}
					}				
			 	}
			 	else if (!theta_i[p] && theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_j*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						if (flags[U_i[j]*(GCA->LUT_size) + q])
						{
							pp = (((q & mask_ii) >> GCA->log2s*ii) << GCA->log2s);
							pp |= ((q & mask_j) >> r); 
							flags[U_i[j]*(GCA->LUT_size) + q] = (p == pp) ? 0 : flags[U_i[j]*(GCA->LUT_size) + q];
						}
					}
			 	}
			}
		}/*end for link i,j*/
	}/*end for cell_i*/
	for (i=0;i<startcell+1;i++)
	{

		U_i = GetNeighbourhood(GCA,i);
		for (j=0;j<(GCA->params->k-1);j++)
		{
			U_j = GetNeighbourhood(GCA,U_i[j]);
			/*get the neighbour number for each cell in the other neighbourhood*/
			ii = 0;
			while (U_j[ii] != i) ii++;
			jj = 0;
			while (U_i[jj] != U_i[j]) jj++;
			ii = (ii >= (GCA->params->k-1)/2) ? ii+1 : ii;
			jj = (jj >= (GCA->params->k-1)/2) ? jj+1 : jj;
			/*create masks for selecting respective Pre-neighbourhood states*/
			mask_i = (state_mask << r); 
			mask_jj = (state_mask << GCA->log2s*jj);
			mask_j = (state_mask << r); 
			mask_ii = (state_mask << GCA->log2s*ii);
				
			/*create list for sub-configs for match test*/
			memset((void*)theta_i,0,(GCA->params->s)*(GCA->params->s)*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/
				if( flags[i*(GCA->LUT_size) + q])
				{
					p = (((q & mask_i) >> r) << GCA->log2s);
					p |= ((q & mask_jj) >> GCA->log2s*jj); 
					/*only set if this still a possible pre-neighbourhood*/
					theta_i[p] = 1;
				}
			}
				
			memset((void*)theta_j,0,(GCA->params->s)*(GCA->params->s)*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/						
				if (flags[U_i[j]*(GCA->LUT_size) + q])
				{
					p = (((q & mask_ii) >> GCA->log2s*ii) << GCA->log2s);
					p |= ((q & mask_j) >> r); 
					/*only set if this still a possible pre-neighbourhood*/
					theta_j[p] = 1;
				}
			}
				
			/* if theta_i not in theta_j then matching Pre-neighbourhoods 
			 * cannot contribute Pre-Image around cell_i, likewise for 
			 * theta_j not in theta_i 
			 * TODO: could repace theta_i and theta_j with bit strings
			 */
			for (p=0;p<(GCA->params->s)*(GCA->params->s);p++)
			{
			 	if (theta_i[p] && !theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_i*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						if (flags[i*(GCA->LUT_size) + q])
						{
							pp = (((q & mask_i) >> r) << GCA->log2s);
							pp |= ((q & mask_jj) >> GCA->log2s*jj); 
							flags[i*(GCA->LUT_size) + q] = (p == pp) ? 0 : flags[i*(GCA->LUT_size) + q];
						}
					}				
			 	}
			 	else if (!theta_i[p] && theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_j*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						if (flags[U_i[j]*(GCA->LUT_size) + q])
						{
							pp = (((q & mask_ii) >> GCA->log2s*ii) << GCA->log2s);
							pp |= ((q & mask_j) >> r); 
							flags[U_i[j]*(GCA->LUT_size) + q] = (p == pp) ? 0 : flags[U_i[j]*(GCA->LUT_size) + q];
						}
					}
			 	}
			}
		}/*end for link i,j*/
	}/*end for cell_i*/
}

/* GetFlags(): creates an array F in N x s^k where element F(i,j) = 0 means neighbourhood
 *             i cannot occur for cell j in any pre-image.
 */
unsigned char* GetFlags(GraphCellularAutomaton* GCA)
{
	unsigned char exit,iter,*flags,*tmp_flags;
	state *theta_i,*theta_j;
	unsigned int i,j,k,ii,jj;
	unsigned int sum,prev_sum;
	unsigned char *counts;
	if (!(flags = (unsigned char*)malloc((GCA->params->N)*(GCA->LUT_size)*sizeof(unsigned char))))
	{
		return NULL;
	}
	
	memset((void*)flags,0,(GCA->params->N)*(GCA->LUT_size)*sizeof(unsigned char));
	if (!(tmp_flags = (unsigned char*)malloc((GCA->params->N)*(GCA->LUT_size)*sizeof(unsigned char))))
	{
		return NULL;
	}
	
	if (!(theta_i = (state*)malloc((GCA->params->s)*(GCA->params->s)*sizeof(state))))
	{
		return NULL;
	}
	
	if (!(theta_j = (state*)malloc((GCA->params->s)*(GCA->params->s)*sizeof(state))))
	{
		return NULL;
	}
	
	if (!(counts = (unsigned char*)malloc((GCA->params->N)*sizeof(unsigned char))))
	{
		return NULL;
	}

	/*TODO: Implement my algorithm for computing all pre-images at once*/
	/* So far I can't see why it would not work, could be a paper in it?*/
	
	/*initial flags - set to 1 if a possible pre-neighbourhood else 0*/
	for (i=0;i<GCA->params->N;i++)
	{
		state s = GetCellStatePacked(GCA,i,0);
		for (j=0;j<GCA->LUT_size;j++)
		{
			flags[i*(GCA->LUT_size)+j] = (GCA->ruleLUT[j] == s) ? 1 : 0;
		}
	}
	prev_sum = 0;
	for (i=0;i<GCA->params->N;i++)
	{
		sum += flags[i];
	}
	
	exit = 0;
	while (!exit)
	{
		iter = 1;
		while (iter)
		{
			RevAlgIter(GCA,flags,theta_i,theta_j,0);
			/*test for convergence*/
			sum = 0;
			for (i=0;i<(GCA->params->N)*(GCA->LUT_size);i++)
			{
				sum += flags[i];
			}
		
			if (sum == 0)
			{
				iter = 0;
				exit = 1;
			}
			/*loop invarient: 0 <= sum <= prev_sum*/
			if (sum == prev_sum)
			{
	
				iter = 0;
			}
			else
			{
				prev_sum = sum;
			}
		}
		/*test of Garden-of-eden states*/
		/*TODO: although colzeros => GOE*/
		/*      !(GOE => colzeros)*/
		/*DAMMIT! gotta rethink here...*/
		for (i=0;i<(GCA->params->N);i++)
		{
			counts[i] = 0;
			for (j=0;j<(GCA->LUT_size);j++)
			{
				counts[i] += flags[i*(GCA->LUT_size)+j];
			}
		}

		for (i=0;i<(GCA->params->N);i++)
		{
			if (counts[i] == 0)
			{
				/*then set the entire flags array to zero*/
				/*this is definitely a GOE*/
				for (j=0;j<(GCA->params->N)*(GCA->LUT_size);j++)
				{
					flags[j] = 0;
				}
				exit = 1;	
			}
			
		}
		/*if we are still uncertain we must further analyse*/
		if(!exit)
		{
			unsigned char invalid;
			invalid = 0;
			memcpy(tmp_flags,flags,(GCA->params->N)*(GCA->LUT_size)*sizeof(unsigned char));
				
			/*for each possible n-hood, run a single test iteration*/
			for (i=0;i<(GCA->params->N);i++)
			{
				for (j=0;j<GCA->LUT_size;j++)
				{
					if( tmp_flags[i*(GCA->LUT_size)+j] != 0)
					{
						/*consider this nhood as fixed*/
						for (k=0;k<j;k++)
						{
							tmp_flags[i*(GCA->LUT_size)+k] = 0;
						}
						for (k=j+1;k<(GCA->LUT_size);k++)
						{
							tmp_flags[i*(GCA->LUT_size)+k] = 0;
						}
						/*run an iteration*/
						RevAlgIter(GCA,tmp_flags,theta_i,theta_j,i);

						/*did it get eliminated*/
						if(tmp_flags[i*(GCA->LUT_size)+j] == 0)
						{
							flags[i*(GCA->LUT_size)+j] = 0;
							invalid = 1;
						}
						else
						{
							memcpy(tmp_flags,flags,(GCA->params->N)*(GCA->LUT_size)*sizeof(unsigned char));
						}
					}

					if (invalid)
					{
						break;
					}
				}

				if(invalid)
				{
					break;
				}
			}
			/*if we get here and we are still valid then we exit*/
			exit = !invalid;
		}
	}
	
	free(theta_i);
	free(theta_j);
	free(tmp_flags);
	free(counts);
	/*flags along with the center state of LUT indexes encodes all possible
	 * pre-images for the current configurations
	 */
	return flags;
};

/* IsGOE(): determines if The current CA configuration is a 
 *          Garden-of-Eden (GOE) configuration.
 * 
 * Paramters:
 *     GCA - graph cellular automaton
 * Returns:
 *     1 if the CA is at a GOE, ottherwise 0
 */
unsigned char IsGOE(GraphCellularAutomaton *GCA)
{
	unsigned char *flags;
	unsigned int sum,i,nbytes;
	nbytes = (GCA->params->N)*(GCA->LUT_size);
	/*Get output from the reverse algorithm pre-processing*/
	flags = GetFlags(GCA);
	sum = 0;
	for (i=0;i<nbytes;i++) sum += flags[i];
	free(flags);
	return sum == 0;
}
/* isValid(): tests if a configuration is a pre-image
 */
unsigned char isValid(GraphCellularAutomaton *GCA,chunk *config,unsigned char *flags)
{
	int i;
	unsigned int nhood;
	for (i=0;i<GCA->params->N;i++)
	{
		nhood = GetNeighbourhood_config_external(GCA,config,i);
		if (!flags[i*GCA->LUT_size + nhood])
		{
			return 0;
		}
	}
	return 1;
}

/* ShannonEntropy(): Computes the Shannon entropy for the CA's spatio-temporal
 *                   pattern. 
 *
 * Parameters:
 *     GCA - the Graph Cellular Automaton
 *     T - number of timesteps to approximate probabilities
 *     pr - storage for probabilities, set to NULL if not required
 *
 * Note:             
 *   The Shannon entropy is defined as
 *        S = 1/N * sum_{i=1}^{N}{sum_{j=0}^{s}{-p_i^j * log_s(p_i^j)}}
 *   Where p_i^j is the probability of cell i having a state j. 
 *   To compute the probabilities the first WSIZE timesteps are omitted, the state counts are
 *   then accumulated accross the number of timesteps.  
 *   See C. Marr et. al.
 */
float ShannonEntropy(GraphCellularAutomaton *GCA, unsigned int T,float *pm,float* logs_pm,float *S_im,unsigned char *TFm,unsigned int *cm)
{
	float *p,*logs_p,*S_i,S,T_inv,logs_inv;
	unsigned char *TF;
	unsigned int *count,t,i,j,N,WSIZE,s;
	
	N = GCA->params->N;
	s = (unsigned int)(GCA->params->s);
	WSIZE = GCA->params->WSIZE;
	
	/*for pr not NULL we assume that pr is correct size*/	
	if(pm != NULL)
	{
		p = pm;
	}
	else
	{
		p = (float *)malloc(N*s*sizeof(float));
		if (!p)
		{
			return -1.0;
		}
	}
	memset((void*)p,0,N*s*sizeof(float));
	
	if (cm != NULL)
	{
		count = cm;
	}
	else
	{
		count = (unsigned int *)malloc(N*s*sizeof(unsigned int));
		if (!count)
		{
			return -1.0;
		}
	}
	memset((void*)count,0,N*s*sizeof(unsigned int));
	
	if (TFm != NULL)
	{
		TF = TFm;
	}
	else
	{
		TF = (unsigned char*)malloc(N*sizeof(unsigned int));
		if (!TF)
		{
			return -1.0;
		}
	}

	if (logs_pm != NULL)
	{
		logs_p = logs_pm;
	}
	else
	{
		logs_p = (float *)malloc(N*s*sizeof(float));
		if (!logs_p)
		{
			return -1.0;
		}
	}
	if (S_im != NULL)
	{
		S_i = S_im;
	}
	else
	{
		S_i = (float *)malloc(N*sizeof(float));
		if (!S_i)
		{
			return -1.0;
		}
	}
	memset((void*)S_i,0,N*sizeof(float));
	
	/*run lead in period*/
	CASimTSteps(GCA,WSIZE);
	
	/*now accumulate probabilities*/
	for (t=0;t<T;t++)
	{
		CANextStep(GCA);
		/*accumulate probabilities*/
		for (j=0;j<s;j++)
		{
			for (i=0;i<N;i++)
			{
				register state s_i = GetCellStatePacked(GCA,i,0);
				TF[i] = (s_i == j);
			}
			for (i=0;i<N;i++)
			{
				count[j*N + i] += TF[i];
			}
		}
	}
	
	T_inv = 1.0/((float)T);
	for (i=0;i<N*s;i++)
	{
		p[i] = T_inv*((float)count[i]);
	}
	
	/*compute log base s for each p_i^j*/
	logs_inv = 1.0/log(s);	
	for (i=0;i<N*s;i++)
	{
		logs_p[i] = logs_inv;
	}
	for (i=0;i<N*s;i++)
	{
		logs_p[i] *= log(p[i]);
	}
	
	/*compute each individual entropy*/
	for (j=0;j<s;j++)
	{
		for (i=0;i<N;i++)
		{
			if (p[j*N+i] != 0.0)
			{
				S_i[i] -= p[j*N+i]*logs_p[j*N+i];
			}
		}
	}
	
	/*accumulate final entropy*/
	S = 0.0;
	for (i=0;i<N;i++)
	{
		S += S_i[i];
	}
	/*average over all cells*/
	S /= (float)N;
	
	/*clean up*/
	if (TFm == NULL)
	{
		free(TF);
	}
	if (logs_pm == NULL)
	{
		free(logs_p);
	}
	if (S_im == NULL)
	{
		free(S_i);
	}
	if (cm == NULL)
	{
		free(count);
	}
	/*only free p if pr was not passed as a parameter*/
	if (pm == NULL)
	{
		free(p);
	}
	return S;
}

/* WordEntropy(): Computes the Word Entropy of the CA's spatio-temporal pattern
 *
 * Parameters:
 *     GCA - The Graph Cellular Automaton
 *     T - the number of timesteps to approximate probabilities
 *     pr - storage for probabilities, set to NULL if not required
 *
 * Note:             
 *   The Word entropy is defined as
 *        W = 1/N * sum_{i=1}^{N}{sum_{l=0}^{T}{-p_i^l * log_s(p_i^l)}}
 *   Where p_i^l is the probability of cell i having a constant word of length l. 
 *   To compute the probabilities the first WSIZE timesteps are omitted, 
 *   the state counts are then accumulated accross the number of timesteps.
 *   See C. Marr et. al.
 */
float WordEntropy(GraphCellularAutomaton *GCA,unsigned int T, float *pm,float* logs_pm,float *W_im,unsigned char *TFm, unsigned int *cm,unsigned int *wlm)
{
	float *p,*logs_p,*W_i,W,T_inv,logs_inv;
	unsigned char *TF;
	unsigned int *count,*word_lengths,i,l,N,WSIZE,s;
	
	N = GCA->params->N;
	WSIZE = GCA->params->WSIZE;
	s = GCA->params->s;
	
	if (pm != NULL)
	{
		p = pm;
	}
	else
	{
		p = (float *)malloc(N*T*sizeof(float));
		if (!p)
		{
			return -1.0;
		}
	}
	memset((void*)p,0,N*T*sizeof(float));

	if (logs_pm != NULL)
	{
		logs_p = logs_pm;
	}
	else
	{
		logs_p = (float *)malloc(N*T*sizeof(float));
		if (!logs_p)
		{
			return -1.0;
		}
	}

	if (W_im != NULL)
	{
		W_i = W_im;
	}
	else
	{
		W_i = (float *)malloc(N*sizeof(float));
		if (!W_i)
		{
			return -1.0;
		}
	}
	memset((void*)W_i,0,N*sizeof(float));

	if (TFm != NULL)
	{
		TF = TFm;
	}
	else
	{
		TF = (unsigned char*)malloc(N*sizeof(unsigned char));
		if (!TF)
		{
			return -1.0;
		}
	}
	memset((void*)TF,0,N*sizeof(unsigned char));
	
	if (wlm != NULL)
	{
		word_lengths = wlm;
	}
	else
	{
		word_lengths = (unsigned int*)malloc(N*sizeof(unsigned int));
		if (!word_lengths)
		{	
			return -1.0;
		}
	}

	if (cm != NULL)
	{
		count = cm;
	}
	else
	{
		count = (unsigned int *)malloc(N*T*sizeof(unsigned int));
		if (!count)
		{
			return -1.0;
		}
	}
	memset((void*)count,0,N*T*sizeof(unsigned int));
	/*initialise word lengths at 1*/
	for (i=0;i<N;i++)
	{
		word_lengths[i] = 1;
	}
	
	/*compute lead in*/
	CASimTSteps(GCA,WSIZE);
	
	/*approximate probabilities*/
	for (l=1;l<T;l++)
	{
		CANextStep(GCA);
		/*check which cells have changed state*/
		for (i=0;i<N;i++)
		{
			TF[i] = (GetCellStatePacked(GCA,i,0) == GetCellStatePacked(GCA,i,1));
		}

		/*increment word lengths for those that did not change*/
		/*set to 1 those that did*/
		for (i=0;i<N;i++)
		{
			word_lengths[i] *= TF[i];
		}
		for (i=0;i<N;i++)
		{
			word_lengths[i]++;
		}
		
		/*increment counts for words that did change*/
		for (i=0;i<N;i++)
		{
			count[word_lengths[i]*N+ i] += TF[i]; 
		}		
	}

	T_inv = 1.0/((float)T);
	for (i=0;i<N*T;i++)
	{
		p[i] = T_inv*((float)count[i]);
	}
	
	/*compute log base s for each p_i^j*/
	logs_inv = 1.0/log(T);	
	for (i=0;i<N*T;i++)
	{
		logs_p[i] = logs_inv;
	}
	for (i=0;i<N*T;i++)
	{
		logs_p[i] *= log(p[i]);
	}
	
	/*compute each individual entropy*/
	for (l=0;l<T;l++)
	{
		for (i=0;i<N;i++)
		{
			if (p[l*N+i] != 0.0)
			{
				W_i[i] -= p[l*N+i]*logs_p[l*N+i];
			}
		}
	}
	
	/*accumulate final entropy*/
	W = 0.0;
	for (i=0;i<N;i++)
	{
		W += W_i[i];
	}
	/*average over all cells*/
	W /= (float)N;
	
	/*clean up*/
	if (TFm == NULL)
	{
		free(TF);
	}
	if (logs_pm == NULL)
	{
		free(logs_p);
	}
	if (W_im == NULL)
	{
		free(W_i);
	}
	if (cm == NULL)
	{
		free(count);
	}
	if (wlm == NULL)
	{
		free(word_lengths);
	}
	/*only free p if pr was not passed as a parameter*/
	if (pm == NULL)
	{
		free(p);
	}

	return W;
}

/* SumCAImages(): Accumulates the counts of each state j per cell i, over all
 *                images of the given pre-Images.
 *
 * Parameters:
 *      GCA - graph cellular automaton
 *      counts - array of length N x s to store cell counts
 *      preImages - an array of pre-Images, or an interval of pre-images
 *      n - number of pre-images, if n == 0 then pre-images contains a lower
 *          and upper range.
 */
unsigned int *SumCAImages(GraphCellularAutomaton *GCA,unsigned int *counts,chunk *preImages,unsigned int n)
{
	unsigned int N,s,i,j,t,p,WSIZE;
	N = GCA->params->N;
	s = (unsigned int)GCA->params->s;
	WSIZE = GCA->params->WSIZE;	
	for (i=0;i<N*s;i++)
	{
		counts[i] = 0;
	}

	/*in this case we have a set list of ICs to test*/
	if (n != 0)
	{
		/*for each pre-image*/
		for (p=0;p<n;p++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,preImages+p*(GCA->size),0);
			/*Run until the time equals the window size*/
			GCA->t = 0;
			t = 0;
			while ((t<WSIZE) && !IsAttCyc(GCA))
			{
				/*step in time*/
				CANextStep(GCA);
				/*accumulate counts*/
				for (j=0;j<s;j++)
				{
					for (i=0;i<N;i++)
					{
						register state s_i = GetCellStatePacked(GCA,i,0);
						counts[j*N + i] += (s_i == j);
					}
				}
				t = GCA->t;
			}
		}
	}
	else
	{
		/*only a single chunk size is legal here for compute reasons*/
		if (GCA->size != 1)
		{
			return NULL;
		}
		/*need to count by two's to avoid rotational symmetry*/	
		for (p=preImages[0];p<preImages[1];p+=2)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,&p,0);
			GCA->t = 0;
			/*only run the simulation if p is a GOE*/
			if (IsGOE(GCA))
			{
				/*Run until the time equals the window size*/
				t = 0;
				while ((t<WSIZE) && !IsAttCyc(GCA))
				{
					/*step in time*/
					CANextStep(GCA);
					/*accumulate counts*/
					for (j=0;j<s;j++)
					{
						for (i=0;i<N;i++)
						{
							register state s_i = GetCellStatePacked(GCA,i,0);
							counts[j*N + i] += (s_i == j);
						}
					}
					t = GCA->t;
				}
			}
		}
	}

	return counts;
}

#ifndef NO_THREADS
/* ComputeExactProbs_worker(): thread entry function for ComputeExactProbs()
 *
 * Parameters:
 *     params - void pointer can be cast to a ThreadData pointer containing a
 *              Graph CA, a range of pre-Images and an array to store counts
 */
void *ComputeExactProbs_worker(void * params)
{
	ThreadData *tData;
	GraphCellularAutomaton *GCA;
	unsigned int *counts;
	chunk *range;
	
	/*unpack data*/
	tData = (ThreadData*) params;
	GCA = (GraphCellularAutomaton *) tData->param1;
	range = (chunk *) tData->param2;
	counts = (unsigned int *)tData->param3;
	
	SumCAImages(GCA,counts,range,0);
	pthread_exit((void*)((unsigned long long int)tData->ID));
}
#endif

/* ComputeExactProbs(): Computes the exact probability of a particular configuration
 *                      occuring after t = 0.
 * Parameters:
 *		GCA - graph cellular automaton
 *
 * Note:
 *    by default this function is threaded using pthreads, if this is desired to 
 *    be disabled then compile with the -DNO_THREADS compiler flag
 */
float *ComputeExactProbs(GraphCellularAutomaton *GCA)
{
	unsigned int i,j,N,s;
	float *probs;
	
#ifndef NO_THREADS	
	GraphCellularAutomaton *GCAs[NUM_THREADS];
	unsigned int *counts[NUM_THREADS];
	chunk *ranges[NUM_THREADS];
	unsigned int work_per_thread,rem;
	/*pthread stuff*/
	pthread_t workers[NUM_THREADS];
	pthread_attr_t attr;
	void *status;
	int r_code;
	ThreadData *tData[NUM_THREADS];
	
#else
	unsigned int *counts;
	chunk range[2];
#endif

	N = GCA->params->N;
	s = (unsigned int)GCA->params->s;
	
#ifndef NO_THREADS	
	/*allocate memory*/
	for (i=0;i<NUM_THREADS;i++)
	{
		tData[i] = (ThreadData*)malloc(sizeof(ThreadData));
		if (!tData[i])
		{
			return NULL;
		}
		
		GCAs[i] = CopyGCA(GCA);
		if (!GCAs[i])
		{
			return NULL;
		}
		
		counts[i] = (unsigned int*)malloc(N*s*sizeof(unsigned int)); 
		if (!counts[i])
		{
			return NULL;
		}
		
		ranges[i] = (chunk *)malloc(2*sizeof(chunk));
		if (!ranges[i])
		{
			return NULL;
		}
	}
#else
	counts = (unsigned int*)malloc(N*s*sizeof(unsigned int)); 
	if (!counts)
	{
		return NULL;
	}
#endif

	probs = (float*)malloc(N*s*sizeof(float)); 
	if (!probs)
	{
		return NULL;
	}
	memset((void *)probs,0,N*s*sizeof(float));
	
#ifndef NO_THREADS	
	/* Initialize and set thread detached attribute */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	work_per_thread = ((0x1 << ((GCA->params->N)*(GCA->log2s))) - 1)/NUM_THREADS;
	rem = ((0x1 << ((GCA->params->N)*(GCA->log2s))) - 1)%NUM_THREADS;
	ranges[0][0] = 0;
	/*first thread takes the slack*/
	ranges[0][1] = work_per_thread + rem;
	for (i=1;i<NUM_THREADS;i++)
	{
		ranges[i][0] = ranges[i-1][1];
		ranges[i][1] = ranges[i][0] + work_per_thread;
	}
	
	
	/*fire up threads...start working them cores*/
	for (i=0;i<NUM_THREADS;i++)
	{	
		/*pack data for each thread*/
		tData[i]->ID = i;		
		tData[i]->param1 = (void *)GCAs[i];
		tData[i]->param2 = (void *)ranges[i];
		tData[i]->param3 = (void *)counts[i];
		
		/*create thread*/
		r_code = pthread_create(&workers[i],&attr,ComputeExactProbs_worker,(void *)tData[i]);
		if(r_code)
		{
			fprintf(stderr,"Thread %d failed to start Return code [%d]\n",i,r_code);
			exit(1);
		}
	}
	
	/* Free attribute and wait for the other threads */
	pthread_attr_destroy(&attr);
	for(i=0;i<NUM_THREADS;i++) 
	{
		r_code = pthread_join(workers[i], &status);
		if(r_code)
		{
			fprintf(stderr,"Thread %d failed to join. Return code [%d]\n",i,r_code);
			exit(1);
		}
	}

	/*accumulate all results*/
	for (i=0;i<NUM_THREADS;i++)
	{
		for (j=0;j<N*s;j++)
		{
			probs[j] += (float)counts[i][j];
		}
	}
#else
	range[0] = 0x0;
	range[1] = (0x1 << ((GCA->params->N)*(GCA->log2s))) - 1;
	SumCAImages(GCA,counts,range,0);
	for (j=0;j<N*s;j++)
	{
		probs[j] += (float)counts[j];
	}
#endif
	
	for (j=0;j<N*s;j++)
	{
		probs[j] /= (float)N;
	}
	
#ifndef NO_THREADS	
	/*clean up*/
	for (i=0;i<NUM_THREADS;i++)
	{
		free(tData[i]);
		free(counts[i]);
		free(GCAs[i]);
		free(ranges[i]);
	}
#else
	free(counts);
#endif

	return probs;
}

/* InputEntropy(): Computes the Input Entropy of the CA over T time-steps, each
 *                 averaged over the WSIZE window.
 * Parameters:
 *		GCA - yet again that Graph CA
 *      T - number of time steps to compute input entropy over
 *      mu - address to store mean input entropy
 *      sigma - address to store input entropy variance
 *
 * Returns:
 *    pointer to array of length T, which stores the input entropy values for each
 *    time-step
 *
 * Note: Input entropy is the Shannon entropy of the lookup table fequency histogram
 *       taken over WSIZE timesteps. See A. Weunsche
 */
float* InputEntropy(GraphCellularAutomaton *GCA,unsigned int T,float* mu, float* sigma,unsigned int *Qm, float *logQm,float *IEm)
{
	unsigned int *Q;
	float *logQ;
	unsigned int w;
	unsigned int n;
	unsigned int N;
	unsigned int t,i,j;
	float numLookups;
	float inv_numLookups;
	float inv_logn;
	float *IE;

	/*Q stores the LUT histogram*/
	N = GCA->params->N;
	w = GCA->params->WSIZE;
	n = GCA->LUT_size;
	numLookups = N*w;
	if (Qm != NULL)
	{
		Q = Qm;
	}
	else
	{
		Q = (unsigned int *)malloc(n*sizeof(unsigned int));
		if (!Q)
		{
			return NULL;
		}
	}

	if (logQm != NULL)
	{
		logQ = logQm;
	}
	else
	{
		logQ = (float*)malloc(n*sizeof(unsigned int));
		if (!logQ)
		{
			return NULL;
		}
	}

	if (IEm != NULL)
	{
		IE = IEm;
	}
	else
	{
		IE = (float*)malloc(T*sizeof(float));
		if (!IE)
		{
			return NULL;
		}
	}

	inv_logn = 1.0/log(n);
	inv_numLookups = 1.0/numLookups;
	/*run lead in period*/
	CASimTSteps(GCA,w);
	for (t=0;t<T;t++)
	{
		/*reset the histogram*/
		memset((void*)Q,0,n*sizeof(unsigned int));
		IE[t] = 0.0;
		/*accumulate the LUT histogram*/
		CANextStep(GCA);
		for (j=0;j<w;j++)
		{
			for (i=0;i<N;i++)
			{
				Q[GetNeighbourhood_config(GCA,i,j)]++;
			}
		}

		for (i=0;i<n;i++)
		{
			
			logQ[i] = (Q[i] == 0) ? 0 : log(((float)Q[i])*inv_numLookups);
		}
		for (i=0;i<n;i++)
		{
			logQ[i] *= inv_logn;
		}
		for(i=0;i<n;i++)
		{
			IE[t] -= ((float)Q[i]*inv_numLookups)*logQ[i];
		}
	}

	/*compute the mean inpute entropy*/
	*mu = 0.0;
	for (t=0;t<T;t++)
	{
		*mu += IE[t];
	}
	*mu /= (float)T;
	/*compute the variance input entropy*/
	*sigma = 0.0;
	for (t=0;t<T;t++)
	{
		*sigma += (IE[t] - *mu)*(IE[t] - *mu);
	}
	*sigma /= (float)T;
	*sigma = sqrt(*sigma);
	if (logQm == NULL)
	{
		free(logQ);
	}
	if (Qm == NULL)
	{
		free(Q);
	}
	return IE;
}

/* lambda_param(): Computes Langton's lambda parameter for the given GA rule
 *
 * Parameters:
 *		GCA - ok, I am pretty sure you shouls know what this is by now
 *
 * Returns:
 *     Langton's lambda
 *
 * Note:
 *    Basically ratio of neighbourhoods that yield a non-quiesient state. 
 *    See C. Langton
 */
float lambda_param(GraphCellularAutomaton *GCA)
{
	unsigned char kn;
	unsigned char n;
	unsigned int i;
	float lambda;
	kn = GCA->LUT_size;
	n = 0;
	/*count quiesent state mappings */
	for (i=0;i<kn;i++)
	{
		n += (GCA->ruleLUT[i] == 0);
	}
	lambda = ((float)(kn - n))/((float)kn);
	return lambda;
}

/* Z_param(): Computes Weunsche's Z parameter for the given GA rule
 *
 * Parameters:
 *		GCA - as above
 *
 * Returns:
 *      Weunsche's Z
 *
 * Note:
 *    Defined as the probability that the next cell in the pre-image is deterministic 
 *    See A.  Weunsche
 */
float Z_param(GraphCellularAutomaton *GCA)
{
	float Z_left;
	float Z_right;
	float *R;
	unsigned char k;
	float Rprod;
	unsigned int i,j,p,q;
	unsigned int n;
	k = GCA->params->k;
	R = (float*)malloc(k*sizeof(float));
	/*compute R's from right-to-left*/
	/*R_k-i = n_(k-i)/(2^k) where n_k-i = number of deterministic 2^(i+1)-tuples*/
	for (i=0;i<k;i++)
	{
		/*compute n_(k-i)*/
		n = 0;
		for (p=0;p<GCA->LUT_size;p++)
		{
			for (q=p+1;q<GCA->LUT_size;q++)
			{
				/*a0,a1,...,1,*,* -> T & a0,a1,...,0,*,* -> not T*/
				n += (((q >> i)^(p >> i)) == 1) 
					&& (GCA->ruleLUT[p] != GCA->ruleLUT[q]);
			}
		}
		R[k-1-i] = ((float)n)/((float)(GCA->LUT_size));
	}
	/*compute R_k U R_k-1 U R_k-2 U ...*/
	Z_left = R[k-1];
	for (i=0;i<k-1;i++)
	{
		Rprod = 1.0;
		for (j=k-i;j<k;j++)
		{
			Rprod *= (1.0 - R[j]);
		}
		Z_left += R[k-1-i]*Rprod;
	}
	/*compute R's from left-tp-right*/
	for (i=0;i<k;i++)
	{
		/*compute n_(k-i)*/
		n = 0;
		for (p=0;p<GCA->LUT_size;p++)
		{
			for (q=p+1;q<GCA->LUT_size;q++)
			{
				/* *,*,1,...,ak-1,ak -> T & *,*,0,...,0,ak-1,ak -> not T*/
				n += ((((q << i)^(p << i)) & (~((1<<k)-1))) == (1<<(k-1))) 
					&& (GCA->ruleLUT[p] != GCA->ruleLUT[q]);
			}
		}
		R[k-1-i] = ((float)n)/((float)(GCA->LUT_size));
	}
	/*compute R_k U R_k-1 U R_k-2 U ...*/
	Z_right = R[k-1];
	for (i=0;i<k-1;i++)
	{
		Rprod = 1.0;
		for (j=k-i;j<k;j++)
		{
			Rprod *= (1.0 - R[j]);
		}
		Z_right += R[k-1-i]*Rprod;
	}

	return (Z_left > Z_right) ? Z_left : Z_right;
}

/* G_density(): Counts the number of Garden-of-Eden Configurations. That is, those
 *              configurations which have no pre-image, and can only exist as initial
 *              conditions.
 * Parameters:
 *		GCA - go figure 
 */
unsigned int G_density(GraphCellularAutomaton *GCA,chunk* ics,unsigned int n)
{
	unsigned int G;
	chunk i;
	G = 0;
	if (n != 0)
	{
		for (i=0;i<n;i++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,ics+i*(GCA->size),EXPLICIT_IC_TYPE);
			/*Garden-of-Eden test*/
			G += IsGOE(GCA);
		}
	}
	else
	{
		for (i=ics[0];i<ics[1];i++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,&i,EXPLICIT_IC_TYPE);
			/*Garden-of-Eden test*/
			G += IsGOE(GCA);
		}
	}
	return G;
}

