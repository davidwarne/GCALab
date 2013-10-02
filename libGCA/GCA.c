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
 * Last Modified: 18/01/2013
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
 *       v 0.15 (19/11/2012) - i. minor bug fix in Create GCA for totalistic and thresh-hold
 *                                based rule codes
 *       v 0.16 (16/01/2013) - i. added windowsize parameter for CreateECA()
 *       v 0.17 (18/01/2013) - i. Added AttLength() and TransLength() functions
 *       v 0.18 (25/02/2013) - i. Added neighbourhood type option, either Von Neumann or Moore
 *                                neighbourhoods are supported
 *                             ii. Created life based rule type (from Carter Bays' paper)
 *       v 0.19 (01/03/2013) - i. Fixed bug in life rule which was incorrectly extracting 
 *                                environment and fertiliy thresholds.
 *                             ii. added population desity function.
 *                             iii. fixed bug in moore neighbourhood construction (same cells
 *                                  were being included multi[le times)
 *                             iv. fix bug in GetNeighbourhood_config when handling variable 
 *                                 neighbourhood sizes.
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

/** 
 * @brief Creates a topology array from a mesh. 
 * @details If mesh is NULL then a regular 1-dimensional genus-1 topology of \a N cells 
 * is created.
 * 
 * @param m mesh structure to generate the topology from, can be NULL.
 * @param N pointer to number of cells.
 * @param k  pointer to neighbourhood size.
 * 
 * @returns An array of unsigned int of length <em>m->vList->numVerts*k</em> or <em>N*k</em>.
 * @retval NULL Indicates that the Topology could not be constructed.
 */
unsigned int * GenerateTopology(unsigned int *N, unsigned char *k,unsigned char nh_type, mesh *m)
{
	unsigned int *graph;
	unsigned int i,j,ii,jj,r,*e;
	unsigned int N_tmp,k_tmp;
	unsigned int num_nb;
	int *f1,*f2;
	unsigned char tf;
	if (m == NULL)
	{
		N_tmp = (N == NULL) ? DEFAULT_NUM_CELLS : *N;
		k_tmp = (k == NULL) ? DEFAULT_NEIGHBOURHOOD_SIZE : *k;
		
		graph = (unsigned int *)malloc(N_tmp*(k_tmp-1)*sizeof(unsigned int));
		if (!graph)
		{
			return NULL;
		}
		memset((void*)graph,0xFF,N_tmp*(k_tmp-1)*sizeof(unsigned int));
		
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
		switch (nh_type)
		{
			/*only cells thats share an edge are neighbours*/
			case VON_NEUMANN_NEIGHBOURHOOD_TYPE:
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
								tf |= (((f1[j] == f2[jj]) && (f1[j+1] == f2[jj+1])) || ((f1[j] == f2[jj+1]) && (f1[j+1] == f2[jj])));
							}
							tf |=  (((f1[j] == f2[k_tmp-2]) && (f1[j+1] == f2[0])) || ((f1[j] == f2[0]) && (f1[j+1] == f2[k_tmp-2])));

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
							tf |= (((f1[k_tmp-2] == f2[jj]) && (f1[0] == f2[jj+1])) || ((f1[k_tmp-2] == f2[jj+1]) && (f1[0] == f2[jj])));
						}

						tf |= (((f1[k_tmp-2] == f2[k_tmp-2]) && (f1[0] == f2[0])) || ((f1[k_tmp-2] == f2[0]) && (f1[0] == f2[k_tmp-2])));
						if (tf && (ii != i))
							break;
					}	
					graph[i*(k_tmp-1) + k_tmp-2] = ii;
				}
				break;
			/*cells that share vertices are neighbours*/
			case MOORE_NEIGHBOURHOOD_TYPE:
				N_tmp = m->fList->numFaces;
				k_tmp = m->fList->maxVerts*4+1;
				graph = (unsigned int *)malloc(N_tmp*(k_tmp-1)*sizeof(unsigned int));
				if (!graph)
				{
					return NULL;
				}
				memset((void*)graph,0xFF,N_tmp*(k_tmp-1)*sizeof(unsigned int));

				for (i=0;i<N_tmp;i++)
				{
					f1 = GetFace_ptr(m->fList,i);
					num_nb = 0;
					for (j=0;j<(m->fList->maxVerts);j++)
					{
						for (ii=0;ii<N_tmp;ii++)
						{
							f2 = GetFace_ptr(m->fList,ii);
							tf = 0;
							for (jj=0;jj<(m->fList->maxVerts);jj++)
							{
								tf |= (f1[j] == f2[jj]);
							}

							if (tf && (ii != i))
							{
								for (jj=0;jj<num_nb;jj++)
								{
									if (graph[i*(k_tmp-1)+jj] == ii)
										break;
								}
								if (jj == num_nb)
								{
									graph[i*(k_tmp-1) + num_nb] = ii;
									num_nb++;
								}
							}
						}
					}
				}
				break;
		}
	}
	
	*N = N_tmp;
	*k = k_tmp;
	return graph;
}

/**
 * @brief Creates consistent parameters for a Graph Cellular Automata on a given topology.
 * 
 * @param m The Geometry of the domain (which is used to derive topology).
 * @param rule_type The type of rule that the rule code represents.
 * @param rule The rule code.
 * @param ws The window-size (i.e., the number of timesteps to store).
 *
 * @returns A \a CellularAutomatonParameters structure.
 * @retval NULL Failed to create the parameter structure.
 */
CellularAutomatonParameters *CreateCAParams(unsigned char nh_type,mesh *m,state s,unsigned char rule_type, unsigned int rule, unsigned int ws)
{
	CellularAutomatonParameters *params;

	params = (CellularAutomatonParameters *)malloc(sizeof(CellularAutomatonParameters));
	if (!params)
	{
		return NULL;
	}

	params->graph = GenerateTopology(&(params->N),&(params->k),nh_type,m);
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

/**
 * @brief Creates a Graph representation of an Elementary Cellular Automaton (ECA).
 *
 * @param N Number of cells.
 * @param k Neighbour size (i.e. 2*r +1 ).
 * @param rule Wolfram's ECA rule notation.
 *
 * @returns A ECA ready for simulation.
 * @retval NULL Failed to created the ECA. 
 */
GraphCellularAutomaton *CreateECA(unsigned int N,unsigned int k,unsigned int rule,unsigned int ws)
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
	params->WSIZE = (ws == 0) ? DEFAULT_WINDOW_SIZE : ws;
	params->rule = rule;
	params->rule_type = CODE_RULE_TYPE;
	params->graph = GenerateTopology(&N_tmp,&k_tmp,DEFAULT_NEIGHBOURHOOD_TYPE,NULL);
	
	ECA = CreateGCA(params);
	
	return ECA;
}

/**
 * @brief Creates a Graph Cellular Automaton (GCA).
 *
 * @param params Parameters defining CA dynamics and topology.
 *
 * @returns A GCA ready for simulation.
 * @retval NULL Failed to create the GCA.
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
			/* here the rule is composed of 4 states-- three states <,=,>, the state of
			 * interest, then a 8-bit integer representing the threshold
			 */ 
			ref = (state)((GCA->params->rule >> 3*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
			thresh = (unsigned char)((GCA->params->rule >> 4*(GCA->log2s)) & 0xFF); 
			for (i=0;i<GCA->LUT_size;i++)
			{
				register unsigned char sum;
				register state s_i;
				register unsigned char ii;
				sum = 0;
				for (j=0;j<GCA->params->k;j++)
				{
					s_i = (state)((i >> j*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
					sum += (s_i == ref);
				}
				
				if (sum < thresh)
				{
					ii = 0;
				}
				else if (sum == thresh)
				{
					ii = 1;
				}
				else
				{
					ii = 2;
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
					s_i = (state)((i >> j*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
					sum += (s_i == ref);
				}
				for (j=0;j<=GCA->params->k;j++)
				{
					if (sum == (state)j)
					{
						ii = (unsigned char)j;
						break;
					}
				}
				
				GCA->ruleLUT[i] = (state)((GCA->params->rule >> (ii+1)*(GCA->log2s)) & ((0x1 << (GCA->log2s)) - 1));
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
		case LIFE_RULE_TYPE: /*Based on Carter Bay's paper*/
		{
			unsigned int E_l,E_h,F_l,F_h;
			unsigned int pop,q,r;
			unsigned int C;
			GCA->LUT_size = pow(params->s,params->k);
			GCA->ruleLUT = (state*)malloc((GCA->LUT_size)*sizeof(state));

			if (!(GCA->LUT_size))
			{
				return NULL;
			}

			/*for a life rule we assume only binary CA*/
			if (GCA->params->s != 2)
			{
				return NULL;
			}

			r = params->rule % 10;
			q = params->rule / 10;
			F_h = r;
			r = q % 10;
			q = q / 10;
			F_l = r;
			r = q % 10;
			q = q / 10;
			E_h = r;
			r = q % 10;
			q = q / 10;
			E_l = r;
			/*life rule has for parts E_l E_h F_l F_h*/
			for (i=0;i<GCA->LUT_size;i++)
			{
				C = (i >> ((GCA->params->k-1)/2)) & 0x1;
				pop = 0;
				for (j=0;j<GCA->params->k;j++)
				{
					pop += (i >> j) & 0x1;
				}
				pop -= C;
				GCA->ruleLUT[i] = (((C==1) && ((pop >= E_l) && (pop <= E_h))) || ((C==0) && ((pop >= F_l) && (pop <= F_h))));
			}
		}
			break;
	}
	
	return GCA;
}

/**
 * @brief Creates copy of the given Graph Cellular Automaton.
 *
 * @param GCA The Graph Cellular Automaton to copy.
 *
 * @returns A GCA identical to the given one.
 * @retval NULL Failed to make a copy of \a GCA.
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

/**
 * @brief Set Cellular Automaton intitial configuration. 
 *
 * @details If \a ic == NULL then initial conditions of the given type are generated.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param ic A pointer to a user specified initial condition.
 * @param type The initial condition type to create.
 *
 * @todo Improve the IC types right now very ECA centric.
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

/**
 * @brief Resets Time-evolution for the Cellular Automaton to \a t0 
 * with the last used initial conditions.
 * 
 * @param GCA The Graph Cellular Automaton to reset
 */
void ResetCA(GraphCellularAutomaton *GCA)
{
	/*reset initial conditions*/
	SetCAIC(GCA,GCA->ic,EXPLICIT_IC_TYPE);
	/*reset time to 0*/
	GCA->t = 0;
}

/**
 * @brief Gets the state of the \a ith cell in a configuration.
 *
 * @details The current time step is referenced by setting \a t = 0. Previous time steps 
 * are referenced by \a t > 0 (e.g., \a t == 1 indicates the previous time step)
 *
 * @param GCA A Graph Cellular Automaton.
 * @param i The index of the cell to extract the state.
 * @param t The time step of interest.
 *
 * @returns The state of the \a ith cell of the GCA at the selected time step.
 *
 * @remark It must be the case that 0 <= \a t < <em>GCA->param->WSIZE</em>
 * @warning This function assumes that \a GCA is using a PACKED cell storage type.
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

/**
 * @brief Sets the state of the \a ith cell in the current configuration to the 
 * given symbol.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param i The index of the cell to update the state of.
 * @param s The new state to assign to cell \a i.
 *
 * @warning This function assumes that \a GCA is using a PACKED cell storage type.
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

/**
 * @brief Gets the state of the ith cell in the given configuration.
 *
 * @details This function does not index the current space-time pattern of \a GCA.
 * Instead in only looks at \a config as though it were temporarily the current 
 * configuration.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param config The configuration.
 * @param i The index of the cell to extract the state of.
 *
 * @returns The state of the \a ith cell in \a config.
 *
 * @warning \a config must be in PACKED format.
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

/**
 * @brief Sets the state of the ith cell in the current configuration to the given 
 * symbol.
 *
 * @details This function will not update the actual current configuration of \a GCA.
 * Instead on the data in \a config is modified.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param config The configuration.
 * @param i The index of the cell to update in \a config.
 *
 * @warning \a config must be in PACKED format.
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

/**
 * @brief Gets the list of cells that are in the neighbourhood of cell i.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param i The index of the cell to query the neighbourhood of.
 *
 * @returns An pointer to an array of cell indexes which are in the neighbourhood 
 * of cell \a i.
 *
 * @warning This function actually returns a reference to the portion of \a
 * GCA graph array, modification of the returned array will corrupt the graph
 * itself (this can of coarse be the desired effect depending on the application).
 */
unsigned int* GetNeighbourhood(GraphCellularAutomaton * GCA,unsigned int i)
{
	return GCA->params->graph + i*(GCA->params->k-1);
}

/**
 * @brief Gets configuration of the neighbourhood of the \a ith cell for the given
 * configuration.
 *
 * @details This function ignores \a GCA's time-space pattern and looks only at the 
 * states of \a config.
 *
 * @param GCA A Graph Cellular Automaton.
 * @param config The configuration.
 * @param i The index of the cell to query the neighbourhood configuration.
 *
 * @returns A bit string (unsigned int) encoding of the local neighbourhood 
 * configuration.
 */
unsigned int GetNeighbourhood_config_external(GraphCellularAutomaton * GCA,chunk* config,unsigned int i)
{
	unsigned int *U_i;
	unsigned int nhood,j,k_local;
	U_i  = GCA->params->graph + i*(GCA->params->k-1);
	/*test if this cell has less neighbours*/
	for (j=0;j<(GCA->params->k-1);j++)
	{
		if (U_i[j] == 0xFFFFFFFF)
		{
			break;
		}
	}
	k_local = GCA->params->k-1;
	nhood = 0;
	nhood |= GetCellStatePacked_external(GCA,config,i) << GCA->log2s*((GCA->params->k-1)/2);
		
	for (j=0;j<(GCA->params->k-1)/2;j++)
	{
		nhood |= GetCellStatePacked_external(GCA,config,U_i[j]) << GCA->log2s*j;
		
	}
	for (j=(GCA->params->k-1)/2;j<k_local;j++)
	{
		nhood |= GetCellStatePacked_external(GCA,config,U_i[j]) << GCA->log2s*(j+1);
	}
	return nhood;
}

/**
 * @brief Gets configuration of the neighbourhood of the \a ith cell.
 *
 * @details The current timestep is referenced by setting \a t = 0. Previous time steps
 * are referenced by \a t > 0 (e.g., \a t == 1 indicates the previous time step).
 *
 * @param GCA A Graph Cellular Automaton.
 * @param i The index of the cell to query the neighbourhood configuration.
 * @param t The time step of interest.
 *
 * @returns A bit string (unsigned int) encoding of the local neighbourhood 
 * configuration.
 */
unsigned int GetNeighbourhood_config(GraphCellularAutomaton * GCA,unsigned int i,unsigned int t)
{
	unsigned int * U_i;
	unsigned int nhood,j,k_local;
	U_i  = GCA->params->graph + i*(GCA->params->k-1);
	/*test if this cell has less neighbours*/
	for (j=0;j<(GCA->params->k-1);j++)
	{
		if (U_i[j] == 0xFFFFFFFF)
		{
			break;
		}
	}
	k_local = j;
	nhood = 0;
	nhood |= GetCellStatePacked(GCA,i,t) << GCA->log2s*((GCA->params->k-1)/2);
		
	for (j=0;j<(GCA->params->k-1)/2;j++)
	{
		nhood |= GetCellStatePacked(GCA,U_i[j],t) << GCA->log2s*j;
		
	}
	for (j=(GCA->params->k-1)/2;j<k_local;j++)
	{
		nhood |= GetCellStatePacked(GCA,U_i[j],t) << GCA->log2s*(j+1);
	}
	return nhood;
}

/**
 * @brief Evolves a Cellular Automaton up to time step \a t.
 *
 *
 * @param GCA A pointer to the graph Cellular Automaton.
 * @param t Timestep to terminate at.
 *
 * @remark After calling this function the current configuration of \a
 * GCA will be that of time \a t, even if it had already evolved beyond that.
 */
void CASimTSteps(GraphCellularAutomaton *GCA,unsigned int t)
{
	if (t < GCA->t)
	{
		ResetCA(GCA);
	}
	while(CANextStep(GCA) != t);
}

/**
 * @brief Evolves the Cellular Automaton 1 timestep.
 *
 * @param GCA A pointer to the graph Cellular Automaton.
 *
 * @returns The current time step.
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

/**
 * @brief Simulates the CA until an attractor cycle begins.
 * 
 * @param GCA A Graph Cellular Automaton.
 * @param t Max time step, abort if a cycle is not reached before this step.
 *
 * @returns The first configuration to be repeated, hence the start of the cycle.
 */
chunk* CASimToAttCyc(GraphCellularAutomaton *GCA,unsigned int t)
{
	chunk *cycle;
	unsigned int i,curt;
	cycle = (chunk*)malloc((GCA->size)*sizeof(chunk));
	if (!cycle)
	{
		return NULL;
	}

	while (!IsAttCyc(GCA) && GCA->t < t) 
	{
		CANextStep(GCA);
	}

	for (i=0;i<GCA->size;i++)
	{
		cycle[i] = GCA->config[i];
	}
	return cycle;
}

/**
 * @brief Detects if the CA has entered an attractor cycle.
 *
 * @param GCA Yes, I'll just assume you know what this is...
 * 
 * @returns Returns the length of the cycle if an attractor cycle has begun, 
 * 0 otherwise. 
 */
unsigned char IsAttCyc(GraphCellularAutomaton *GCA)
{
	unsigned int WSIZE,nbytes,t,tn;
	WSIZE = GCA->params->WSIZE;
	nbytes = (GCA->size)*sizeof(chunk);
	if (GCA->t < WSIZE)
	{
		tn = GCA->t;
	}
	else
	{
		tn = WSIZE - 1;
	}
	for (t=tn;t>0;t--)
	{
		/*if we have seen this configuation before, then we have entered an attractor cycle */
		if (!memcmp((void*)(GCA->config),(void*)(GCA->st_pattern[t]),nbytes)){
			/*woah!? dejavu... was it the same cat?*/
			return t;
		}
	}
	return 0;
}

/**
 * @brief Gets all the Pre-Images that can lead to this CAs current configuration.
 * 
 * @param GCA A graph cellular automaton.
 * @param n A reference to store the number of pre-images returned.
 * @param flags The array of consistent local neighbourhood configurations.
 *
 * @returns An array of configurations that are the pre-images of GCA's current 
 * configuration.
 *
 * @retval NULL Configuration is a Garden-of-Eden configuration.
 * 
 * @todo A bit inefficient, and we miss some valid pre-images sometimes
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

/**
 * @brief Implementation of the Neighbourhood Elimination operation. 
 *
 * @details Eliminates currently inconsistent neighbourhoods for each cell.
 * 
 * @param GCA You should definitely know what this by now :).
 * @param flags An array <em>F : F(i,j) = 0</em> => neighbourhood \a i is not possible 
 * for cell \a j in any pre-image.
 * @param theta_i Buffer to store boundary configurations of i to j.
 * @param theta_j Buffer to store boudnary configurations of j to i.
 * @param startcell Last cell to be visited in the algorithm.
 *     
 */
unsigned char NhElim(GraphCellularAutomaton *GCA,unsigned char *flags,state *theta_i,state *theta_j,unsigned int startcell)
{
	unsigned char r;
	unsigned char k2;
	state *W_i,*W_j;
	unsigned int i,j;
	unsigned int sum,prev_sum;
	unsigned int *U_i,*U_j;
	unsigned int state_mask,mask_nh;
    unsigned int theta_size;

	k2 = (GCA->params->k-1)/2;
    r = GCA->log2s*k2;
	
    state_mask = (0x1 << (GCA->log2s)) - 1;
	mask_nh = (state_mask << r); 
	theta_size = (GCA->params->s)*(GCA->params->s);
    /*for every connection*/
	for (i=startcell+1;i<GCA->params->N;i++)
	{
        register unsigned int i_offset;
        i_offset = i*(GCA->LUT_size);
		U_i = GetNeighbourhood(GCA,i);
		for (j=0;j<(GCA->params->k-1);j++)
		{
            register unsigned int ii,jj,q,p,pp;
            register unsigned int mask_ii,mask_jj,log2s_ii,log2s_jj;
			register unsigned int Uij_offset;
            Uij_offset = U_i[j]*(GCA->LUT_size);
			U_j = GetNeighbourhood(GCA,U_i[j]);
			/*get the neighbour number for each cell in the other neighbourhood*/
			ii = 0;
			while (U_j[ii] != i) ii++;
			jj = j;
			ii += (ii >= k2);
			jj += (jj >= k2);
            log2s_ii = GCA->log2s*ii;
            log2s_jj = GCA->log2s*jj;
            /*create masks for selecting respective Pre-neighbourhood states*/
			mask_jj = (state_mask << log2s_jj);
			mask_ii = (state_mask << log2s_ii);
				
			/*create list for sub-configs for match test*/
			memset((void*)theta_i,0,theta_size*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
                /*p = state i : state j*/
				p = (((q & mask_nh) >> r) << GCA->log2s);
				p |= ((q & mask_jj) >> log2s_jj); 
				theta_i[p] |= flags[i_offset + q];

			}
				
			memset((void*)theta_j,0,theta_size*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/						
				p = (((q & mask_ii) >> log2s_ii) << GCA->log2s);
				p |= ((q & mask_nh) >> r); 
				theta_j[p] |= flags[Uij_offset +q];
			}
				
			/* if theta_i not in theta_j then matching Pre-neighbourhoods 
			 * cannot contribute Pre-Image around cell_i, likewise for 
			 * theta_j not in theta_i 
			 */
			for (p=0;p<theta_size;p++)
			{
			 	if (theta_i[p] && !theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_i*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						pp = (((q & mask_nh) >> r) << GCA->log2s);
						pp |= ((q & mask_jj) >> log2s_jj); 
						flags[i_offset + q] &= (p != pp);
					}				
			 	}
			 	else if (!theta_i[p] && theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_j*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						pp = (((q & mask_ii) >> log2s_ii) << GCA->log2s);
						pp |= ((q & mask_nh) >> r); 
						flags[Uij_offset + q] &= (p != pp);
					}
			 	}
			}
		}/*end for link i,j*/
	}/*end for cell_i*/
	for (i=0;i<startcell+1;i++)
	{
        register unsigned int i_offset;
        i_offset = i*(GCA->LUT_size);
		U_i = GetNeighbourhood(GCA,i);
		for (j=0;j<(GCA->params->k-1);j++)
		{
            register unsigned int ii,jj,q,p,pp;
            register unsigned int mask_ii,mask_jj,log2s_ii,log2s_jj;
			register unsigned int Uij_offset;
            Uij_offset = U_i[j]*(GCA->LUT_size);
			U_j = GetNeighbourhood(GCA,U_i[j]);
			/*get the neighbour number for each cell in the other neighbourhood*/
			ii = 0;
			while (U_j[ii] != i) ii++;
			jj = j;
			ii += (ii >= k2);
			jj += (jj >= k2);
            log2s_ii = GCA->log2s*ii;
            log2s_jj = GCA->log2s*jj;
			/*create masks for selecting respective Pre-neighbourhood states*/
			mask_jj = (state_mask << log2s_jj);
			mask_ii = (state_mask << log2s_ii);
				
			/*create list for sub-configs for match test*/
			memset((void*)theta_i,0,theta_size*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/
				p = (((q & mask_nh) >> r) << GCA->log2s);
				p |= ((q & mask_jj) >> log2s_jj); 
				theta_i[p] |= flags[i_offset + q];
			}
				
			memset((void*)theta_j,0,theta_size*sizeof(state));
			for (q=0;q<GCA->LUT_size;q++)
			{
				/*p = state i : state j*/						
				p = (((q & mask_ii) >> log2s_ii) << GCA->log2s);
				p |= ((q & mask_nh) >> r); 
				theta_j[p] |= flags[Uij_offset +q];
			}
				
			/* if theta_i not in theta_j then matching Pre-neighbourhoods 
			 * cannot contribute Pre-Image around cell_i, likewise for 
			 * theta_j not in theta_i 
			 */
			for (p=0;p<theta_size;p++)
			{
			 	if (theta_i[p] && !theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_i*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						pp = (((q & mask_nh) >> r) << GCA->log2s);
						pp |= ((q & mask_jj) >> log2s_jj); 
						flags[i_offset + q] &= (p != pp);
					}				
			 	}
			 	else if (!theta_i[p] && theta_j[p])
			 	{
			 		/*remove invalid pre-neighbourhoods for cell_j*/
			 		for (q=0;q<GCA->LUT_size;q++)
					{
						pp = (((q & mask_ii) >> log2s_ii) << GCA->log2s);
						pp |= ((q & mask_nh) >> r); 
						flags[Uij_offset + q] &= (p != pp);
					}
			 	}
			}
		}/*end for link i,j*/
	}/*end for cell_i*/
}

/**
 * @brief Creates an array <em>F in N x s^k</em> where element <em>F(i,j) = 0</em> 
 * means neighbourhood \a i cannot occur for cell \a j in any pre-image.
 * 
 * @param GCA Well what do you reckon?!!
 *
 * @returns An array such that <em>flags(i,j) = 0 =></em> neighbourhood \a i cannot 
 * occur for cell \a j in any pre-image.
 *
 * @remark This function implements the EDEN-DET algorithm.
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
	
	
	/*initial flags - set to 1 if a possible pre-neighbourhood else 0*/
	for (i=0;i<GCA->params->N;i++)
	{
		state s = GetCellStatePacked(GCA,i,0);
		for (j=0;j<GCA->LUT_size;j++)
		{
			flags[i*(GCA->LUT_size)+j] = (GCA->ruleLUT[j] == s);
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
			NhElim(GCA,flags,theta_i,theta_j,0);
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

        exit = (sum == 0);

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
						NhElim(GCA,tmp_flags,theta_i,theta_j,i);

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

					if (invalid )
					{
						break;
					}
				}

				if(invalid )
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
	/*flags along with the center state of LUT indexes encodes all possible
	 * pre-images for the current configurations
	 */
	return flags;
};

/**
 * @brief Determines if The current CA configuration is a Garden-of-Eden (GOE) 
 * configuration.
 * 
 * @param GCA A graph cellular automaton
 * 
 * @retval 1 if the CA is at a GOE
 * @retval 0 It is very likely CA is not at a GOE.
 *
 * @remark The general Eden problem for graph cellular automata is NP-complete. Thus
 * an exact solution may not exist. It is possible for this function to return 0 but 
 * the CA is actually at a GOE (This occurs is around 15 % of Chaotic CA, but far less 
 * for other dynamical classes).
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
/**
 * @brief Tests if a configuration is a pre-image.
 *
 * @parma GCA A Graph Cellular Automaton.
 * @param config A configuration to test.
 * @param flags An array containing valid neighbourhood configruations.
 *
 * @retval 1 if the CA is at a GOE
 * @retval 0 It is very likely CA is not at a GOE.
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

/**
 * @brief Computes the Shannon entropy for the CA's spatio-temporal pattern. 
 *
 * @param GCA A Graph Cellular Automaton.
 * @param T Number of timesteps to approximate probabilities.
 * @param pm Memory for probablities.
 * @param logs_pm Memory for log of probabilities.
 * @param S_im Memory for Shannon entropy of cells.
 * @param TFm Memory to store state equality test results.
 * @param cm Memory to store counts.
 *
 * @returns The Average Shannon entropy for the CA's Evolution.
 * 
 * @note The Shannon entropy is defined as  <em>S = 1/N * sum_{i=1}^{N}{sum_{j=0}^{s}
 * {-p_i^j * log_s(p_i^j)}}</em> Where <em>p_i^j</em> is the probability of cell \a i having a state \a j. 
 * To compute the probabilities the first WSIZE timesteps are omitted, 
 * the state counts are then accumulated accross the number of timesteps.  
 * See C. Marr et. al.
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

/**
 * @brief Computes the Word Entropy of the CA's spatio-temporal pattern.
 *
 * @param GCA The Graph Cellular Automaton.
 * @param T The number of time steps to approximate probabilities.
 * @param pm Memory for probablities.
 * @param logs_pm Memory for log of probabilities.
 * @param W_im Memory for Word entropy of cells.
 * @param TFm Memory to store state equality test results.
 * @param cm Memory to store counts.
 * @param wlm Memory to store wold lengths.
 *
 * @returns The average Word entropy for the CA's evolution.
 *
 * @note The Word entropy is defined as <em>W = 1/N * sum_{i=1}^{N}{sum_{l=0}^{T}
 * {-p_i^l * log_s(p_i^l)}}</em> Where <em>p_i^l</em> is the probability of cell \a i 
 * having a constant word of length \a l. To compute the probabilities the first 
 * WSIZE timesteps are omitted, the state counts are then accumulated accross the 
 * number of timesteps. See C. Marr et. al.
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
		word_lengths[i] = 0;
	}
	
	/*compute lead in*/
	CASimTSteps(GCA,WSIZE);
	
	/*approximate probabilities*/
	for (l=0;l<T-1;l++)
	{
		CANextStep(GCA);
		/*check which cells have changed state*/
		for (i=0;i<N;i++)
		{
			TF[i] = (GetCellStatePacked(GCA,i,0) == GetCellStatePacked(GCA,i,1));
		}

		for (i=0;i<N;i++)
		{
			if (TF[i])
			{
				/*increment word lengths for those that did not change*/
				word_lengths[i]++;
			}
			else
			{
				/*increment counts for words that did change*/
				count[word_lengths[i]*N+i]++;
				/*set to 1 those that did*/
				word_lengths[i] = 1;
			}

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

/**
 * @brief Accumulates the counts of each state \a j per cell \a i, 
 * over all images of the given pre-Images.
 *
 * @param GCA A Graph cellular automaton.
 * @param counts Array of length \a N x \a s to store cell counts.
 * @param preImages An array of pre-Images, or an interval of pre-images.
 * @param n Number of pre-images, if \a n == 0 then pre-images contains a 
 * lower and upper range.
 *
 * @returns A pointer to an \a N x \a s array, this is equal to the \a counts array.
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
			SetCAIC(GCA,preImages+p*(GCA->size),EXPLICIT_IC_TYPE);
			ResetCA(GCA);
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
			ResetCA(GCA);
			/*only run the simulation if p is a GOE*/
			if (IsGOE(GCA))
			{
				/*Run until the time equals the window size*/
				t = 0;
				while ((t<WSIZE) && !IsAttCyc(GCA))
				{
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

/**
 * @brief Computes the exact probability of a particular configuration occurring after \a t = 0.
 * 
 * @param GCA A graph cellular automaton.
 *
 * @returns An \a N x \a s array of probabilities.
 *
 */
float *ComputeExactProbs(GraphCellularAutomaton *GCA)
{
	unsigned int i,j,N,s;
	float *probs;
	
	unsigned int *counts;
	chunk range[2];

	N = GCA->params->N;
	s = (unsigned int)GCA->params->s;
	
	counts = (unsigned int*)malloc(N*s*sizeof(unsigned int)); 
	if (!counts)
	{
		return NULL;
	}

	probs = (float*)malloc(N*s*sizeof(float)); 
	if (!probs)
	{
		return NULL;
	}
	memset((void *)probs,0,N*s*sizeof(float));
	
	range[0] = 0x0;
	range[1] = (0x1 << ((GCA->params->N)*(GCA->log2s))) - 1;
	SumCAImages(GCA,counts,range,0);
	for (j=0;j<N*s;j++)
	{
		probs[j] += (float)counts[j];
	}
	
	for (j=0;j<N*s;j++)
	{
		probs[j] /= (float)N;
	}
	
	free(counts);

	return probs;
}

/**
 * @brief Computes the Input Entropy of the CA over T time-steps. 
 * 
 * @details The entropy at each time step is averaged over the WSIZE window.
 * @param GCA Yet again that Graph CA.
 * @param T Number of time steps to compute input entropy over.
 * @param mu Address to store mean input entropy.
 * @param sigma Address to store input entropy variance.
 * @param Qm Memory to use for lookup frequencies.
 * @param logQm Memory to use for logarithms of lookup frequencies.
 * @param IEm Memory for the Input Entropy, this pointer is returned as output.
 *
 * @return A pointer to array of length T, which stores the input entropy values for each time-step.
 *
 * @note Input entropy is the Shannon entropy of the lookup table fequency histogram taken over WSIZE timesteps. See A. Weunsche
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
	CASimTSteps(GCA,1000);
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

/**
 * @brief Computes Langton's lambda parameter for the given GA rule.
 *
 * @param GCA OK, I am pretty sure you should know what this is by now.
 *
 * @returns Langton's lambda.
 *
 * @note Basically ratio of neighbourhoods that yield a non-quiesient state. See C. Langton
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

/**
 * @brief Computes Weunsche's Z parameter for the given GA rule.
 *
 * @param GCA As above.
 *
 * @returns Weunsche's Z.
 *
 * @note Defined as the probability that the next cell in the pre-image is deterministic. See A.  Weunsche
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

/**
 * @brief Computes the density of Garden-of-Eden Configurations. 
 *
 * @details Garden-of-Eden configurations are those configurations which have no pre-image, and can only exist as initial conditions.
 * 
 * @param GCA Go figure .
 * @param ics A set of configurations to test, or a range of configurations to test if \a n == 0.
 * @param n If \a ics == \a NULL then this is the number of random samples to use, else it is is the number of configurations in \a ics.
 * 
 * @returns The density of Garden of eden configurations.
 */
float G_density(GraphCellularAutomaton *GCA,chunk* ics,unsigned int n)
{
	unsigned int G;
	chunk i;
	G = 0;
	if (n != 0)
	{
		if (ics != NULL)
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,ics+i*(GCA->size),EXPLICIT_IC_TYPE);
				ResetCA(GCA);
				/*Garden-of-Eden test*/
				G += IsGOE(GCA);
			}
		}
		else
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				ResetCA(GCA);
				/*Garden-of-Eden test*/
				G += IsGOE(GCA);
			}
		}
		return ((float)G)/((float)n);
	}
	else
	{
		for (i=ics[0];i<ics[1];i++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,&i,EXPLICIT_IC_TYPE);
			ResetCA(GCA);
			/*Garden-of-Eden test*/
			G += IsGOE(GCA);
		}
		return ((float)G)/((float)ics[1]-ics[0]);

	}
}

/**
 * @brief Returns the average length of attractor cycles.
 *
 * @param GCA Go figure .
 * @param ics A set of configurations to test, or a range of configurations if \a n == 0.
 * @param n If \a ics == \a NULL then this is the number of random samples to use, else it is is the number of configurations in \a ics.
 * @param t Max time step to simulate before search for an attractor is halted.
 *
 * @returns The average attractor cycle length.
 */
float AttLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned int t)
{
	unsigned int i,numcycles,totaloflengths,length;
	numcycles = 0;
	totaloflengths = 0;
	if ( n != 0)
	{
		if (ics != NULL)
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,ics+i*(GCA->size),EXPLICIT_IC_TYPE);
				ResetCA(GCA);
				/*simulate CA until an Attractor cycle is reached*/
				while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
				{
					CANextStep(GCA);
				}
				totaloflengths += length;
				numcycles += (unsigned int)(length != 0);
			}
		}
		else
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				ResetCA(GCA);
				/*simulate CA until an Attractor cycle is reached*/
				while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
				{
					CANextStep(GCA);
				}
				totaloflengths += length;
				numcycles += (unsigned int)(length != 0);
			}
		}
	}
	else
	{
		for (i=ics[0];i<ics[1];i++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,&i,EXPLICIT_IC_TYPE);
			ResetCA(GCA);
			/*simulate CA until an Attractor cycle is reached*/
			while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
			{
				CANextStep(GCA);
			}
			totaloflengths += length;
			numcycles += (unsigned int)(length != 0);
		}
	}

	return (numcycles > 0) ? ((float)totaloflengths)/((float)numcycles): 0.0;
}

/**
 * @brief Returns the average transient path length.
 *
 * @param GCA Go figure .
 * @param ics A set of configurations to test, or a range of configurations if \a n == 0.
 * @param n If \a ics == \a NULL then this is the number of random samples to use, else it is is the number of configurations in \a ics.
 * @param t Max time step to simulate before search for an attractor is halted.
 * 
 * @returns The average transient path length.
 */
float TransLength(GraphCellularAutomaton *GCA,chunk *ics, unsigned int n,unsigned int t)
{
	unsigned int i,numtrans,totaloflengths,length;
	numtrans = 0;
	totaloflengths = 0;
	if ( n != 0)
	{
		if (ics != NULL)
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,ics+i*(GCA->size),EXPLICIT_IC_TYPE);
				ResetCA(GCA);
				/*simulate CA until an Attractor cycle is reached*/
				while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
				{
					CANextStep(GCA);
				}
				/*test that we did not start within an attractor cycle*/
				if (GCA->t >= length)
				{
					/*we need to subtract the cycle length from the timesteps*/
					totaloflengths += (GCA->t + 1 - length);
					numtrans++;
				}
			}
		}
		else
		{
			for (i=0;i<n;i++)
			{
				/*fix the initial condition*/
				SetCAIC(GCA,NULL,NOISE_IC_TYPE);
				ResetCA(GCA);
				/*simulate CA until an Attractor cycle is reached*/
				while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
				{
					CANextStep(GCA);
				}
				/*test that we did not start within an attractor cycle*/
				if (GCA->t >= length)
				{
					/*we need to subtract the cycle length from the timesteps*/
					totaloflengths += (GCA->t + 1 - length);
					numtrans++;
				}
			}
		}
	}
	else
	{
		for (i=ics[0];i<ics[1];i++)
		{
			/*fix the initial condition*/
			SetCAIC(GCA,&i,EXPLICIT_IC_TYPE);
			ResetCA(GCA);
			/*simulate CA until an Attractor cycle is reached*/
			while(!(length = IsAttCyc(GCA)) && GCA->t < t) 
			{
				CANextStep(GCA);
			}
			/*test that we did not start within an attractor cycle*/
			if (GCA->t >= length)
			{
				/*we need to subtract the cycle length from the timesteps*/
				totaloflengths += (GCA->t + 1 - length);
				numtrans++;
			}
		}
	}

	return (numtrans > 0) ? ((float)totaloflengths)/((float)numtrans) : 0.0;
}

/**
 * @brief Calculates the "live" population density.
 *
 * @details Calculates probability that a cell is non-quiescent over the coarse of the CA's evolution.
 * 
 * @param GCA Go figure...
 * @param ics The initial configuration, if set to NULL then a random initial configuration is used.
 * @param T The number of time steps. 
 * @param dense The array to store densities.
 *
 * @returns An array of probabilities where <em>p[i]</em> is the probability that cell \a i is non-quiescent   .
 */
float* PopDensity(GraphCellularAutomaton *GCA,chunk* ics,unsigned int T, float *dense)
{
	int i,j;
	unsigned int count;
	
	/*test if the user has supplied their own memory refernce*/
	if (!dense)
	{
		dense = (float*)malloc(T*sizeof(float));
		if (!dense)
		{
			return NULL;
		}
	}

	/*test if the user has supplied their own initial conditions*/
	if (ics)
	{
		SetCAIC(GCA,ics,EXPLICIT_IC_TYPE);
		ResetCA(GCA);
	}
	else
	{
		SetCAIC(GCA,NULL,NOISE_IC_TYPE);
		ResetCA(GCA);
	}

	/*simulate and recode population density*/
	for (i=0;i<T;i++)
	{
		count = 0;
		for (j=0;j<(GCA->params->N);j++)
		{
			count += (GetCellStatePacked(GCA,j,0) > 0);
		}
		dense[i] = ((float)count)/((float)GCA->params->N);
		CANextStep(GCA);
	}
	return dense;
}
