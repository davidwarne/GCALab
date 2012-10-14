/* File: test.c
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
 * Description: simple test driver, will eventually be replaced by a proper 
 *              command-line Graph CA simulation tool.
 */

#include <stdio.h>
#include <stdlib.h>
#include "GCA.h"
#define P 50
unsigned char bitcount[16] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
chunk old[1000];

void testbitaccess(int argc, char** argv)
{
	unsigned char D[3];
	unsigned int S,i,lg2,w,sc;
	lg2 = 0;
	S = (unsigned int)atoi(argv[1]);
	i = (unsigned int)atoi(argv[2]);
	D[0] = 0xF6;
	D[1] = 0x75;
	D[2]= 0x43;
	sc = S;
	while(sc >>= 1){
		lg2++;
	}
	printf("%x\n",lg2);
	w = 8/lg2;
	printf("%d %x %d %x\n",i/w,(1<<lg2)-1,i%w,D[i/w]);
	printf("%x\n",(D[i/w] & (((1 << lg2)-1) << (i%w)*lg2)) >> (i%w)*lg2);
}

void testECA(int argc, char **argv)
{
	unsigned int i,j,N,k,rule;
	GraphCellularAutomaton *ECA;
	
	N = (unsigned int)atoi(argv[1]);
	k = (unsigned int)atoi(argv[2]);
	rule = (unsigned int)atoi(argv[3]);
	
	ECA = CreateECA(N,k,rule);
	
	printf("%u %u %u %hhu %hhu %hhu\n",ECA->params->N,ECA->params->WSIZE,ECA->params->rule,ECA->params->s,ECA->params->rule_type,ECA->params->k);
	for (i=0;i<N;i++)
	{
		for (j=0;j<k-1;j++)
		{
			printf(" %u",ECA->params->graph[i*(k-1)+j]);
		}
		printf("\n");
	}
	printf("%hhu %hhu %u %u\n",ECA->log2s,ECA->LUT_size,ECA->size,ECA->t);
	for (i=0;i<ECA->LUT_size;i++)
	{
		printf("%x %hhu\n",i,ECA->ruleLUT[i]);
	}
	ECA->config[0] = 0x00000001;
	for(i=1;i<ECA->size;i++)
		ECA->config[i] = 0x00000000;
	printf("%x\n",ECA->config[0]);
	for (j=0;j<ECA->params->N;j++)
	{
		printf("%u",GetCellStatePacked(ECA,j,0));
	}
	printf("\n");
	for (i=1;i!=100;i++)
	{
		CANextStep(ECA);
		for (j=0;j<ECA->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECA,j,0));
		}
		printf("  %u\n",i);
		
	}
	printf("%x\n",ECA->config[0]);
}

void testRevAlgorithm(int argc,char **argv)
{
	unsigned int i,j,N,k,rule,n;
	unsigned char *flags;
	chunk *preImages;
	GraphCellularAutomaton *ECA;
	
	N = (unsigned int)atoi(argv[1]);
	k = (unsigned int)atoi(argv[2]);
	rule = (unsigned int)atoi(argv[3]);
	
	ECA = CreateECA(N,k,rule);
	
	printf("%u %u %u %hhu %hhu %hhu\n",ECA->params->N,ECA->params->WSIZE,ECA->params->rule,ECA->params->s,ECA->params->rule_type,ECA->params->k);
	ECA->config[0] = 0x1000;
	for (j=0;j<ECA->params->N;j++)
	{
		printf("%u",GetCellStatePacked(ECA,j,0));
	}
	printf("\n");
	flags = GetFlags(ECA);
	
	for (j=0;j<ECA->LUT_size;j++)
	{
		for (i=0;i<ECA->params->N;i++)
		{
			printf("%hhu ",flags[i*ECA->LUT_size+j]);
		}
		printf("\n");
	}
	preImages = CAGetPreImages(ECA,&n,NULL);
	printf("%u\n",n);
	
	for (i=0;i<n;i++)
	{
		for (j=0;j<ECA->params->N;j++)
		{
			printf("%u",GetCellStatePacked_external(ECA,preImages+i*(ECA->size),j));
		}
		printf("\n");
	}

	for (i=0;i<n;i++)
	{
		SetCAIC(ECA,preImages+i*(ECA->size),0);
		printf("test %d\n",i);
		for (j=0;j<ECA->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECA,j,0));
		}
		printf("\n");
		CANextStep(ECA);
		for (j=0;j<ECA->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECA,j,0));
		}
		printf("\n");

	}
	// GetLeastSetBits(ECA,NULL);
	return;
}

void testSE(int argc,char **argv)
{
	unsigned int i,j,N,k,rule;
	GraphCellularAutomaton *ECA;
	
	N = (unsigned int)atoi(argv[1]);
	k = (unsigned int)atoi(argv[2]);
	rule = (unsigned int)atoi(argv[3]);
	
	ECA = CreateECA(N,k,rule);
	
	SetCAIC(ECA,NULL,NOISE_IC_TYPE);
	
	fprintf(stdout,"Shannon Entropy = %f\n",ShannonEntropy(ECA,1000000,NULL));
}
unsigned int popcount(void* bits,int nbytes)
{
	int i=0;
	unsigned int pop;
	char *bits2;
	bits2 = (char *)bits;
	pop = 0;
	for (i=0;i<nbytes;i++){
		printf("\n%u %u %u %u\n",(bits2[i] >> 4) & 0x0F,bits2[i] & 0x0F,bitcount[bits2[i]>>4],bitcount[bits2[i] & 0x0F]);
		pop += bitcount[(bits2[i] >> 4)& 0x0F] + bitcount[bits2[i] & 0x0f];
	}
	return pop;
}
void GetLeastSetBits(GraphCellularAutomaton *GCA,chunk* leastsetbit_config)
{
	chunk *preImages;
	unsigned int n,i,j;
	unsigned int min_set_bits;
	unsigned int min_i;
	preImages = CAGetPreImages(GCA,&n,NULL);
	unsigned int sum = 0;
	min_set_bits = GCA->params->N;
	for (i=0;i<n;i++)
	{
		sum = 0;
		for (j=0;j<GCA->params->N;j++)
		{
			sum +=GetCellStatePacked_external(GCA,preImages + i*GCA->size,j);
		}
		if (sum < min_set_bits)
		{
			min_set_bits = sum;
			min_i = i;
		}
		
	}
	
	if (n > 0)
		//memcpy((void*)leastsetbit_config,(void *)preImages + min_i*GCA->size,GCA->size*CHUNK_SIZE_BITS/8);
		for (i=0;i<GCA->size;i++)
			leastsetbit_config[i] = (preImages + min_i*GCA->size)[i];
	else
		memset((void *)leastsetbit_config,0xff,GCA->size*CHUNK_SIZE_BITS/8);
	return;
}
void GetLeastBits(GraphCellularAutomaton *GCA,unsigned char* flags)
{
	unsigned int i,j,q,p,*U_i,*U_j,ii,jj;
	unsigned int mask_i,mask_j,mask_ii,mask_jj;
		unsigned int state_mask;
	unsigned int rd;
	unsigned char r,accept;
	state *theta_j;
	
	/*create a list to */
	
	if (!(theta_j = (state*)malloc((GCA->params->s)*(GCA->params->s)*sizeof(state))))
	{
		return;
	}
	state_mask = (0x1 << (GCA->log2s)) - 1;
	r = GCA->log2s*(GCA->params->k-1)/2;
	
	for (i=0;i<GCA->params->N;i++)
	{
		
		accept = 0;
		while (!accept)
		{
			unsigned int min = 100000000;
			unsigned int minj,minp;
		/*get max zero neighbourhood*/
		for (j=0;j<GCA->LUT_size;j++)
		{
			if (flags[i*(GCA->LUT_size) + j])
			{
				if (min > bitcount[j])
				{
					min = bitcount[j];
					minj=j;
				}
			}
		}
		if( min==100000000)
		{
			//printf("Error:");
			return;

		}
		/*test that the neighbourhood is possible*/
		U_i = GCA->params->graph + i*(GCA->params->k-1);
		for (j=0;j<(GCA->params->k-1);j++)
		{
			U_j = GCA->params->graph + U_i[j]*(GCA->params->k-1);
				
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
		
			minp = (((minj & mask_i) >> r) << GCA->log2s);
			minp |= ((minj & mask_jj) >> GCA->log2s*jj);
			
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
			
			if (!theta_j[minp])
			{
				
				flags[i*(GCA->LUT_size) + minj] = 0;
			}		
		}
		
		if (flags[i*(GCA->LUT_size) + minj])
		{
			
				
			accept = 1;
			for (j=0;j<minj;j++)
				flags[i*(GCA->LUT_size) + j] = 0;
			for (j=minj+1;j<GCA->LUT_size;j++)
				flags[i*(GCA->LUT_size) + j] = 0;
			}
		}
		
		
	}

	for (i=0;i<GCA->params->N;i++)
	{
		for (j=0;j<GCA->LUT_size;j++)
		{
			if (flags[i*(GCA->LUT_size) + j])
			{
				SetCellStatePacked(GCA,i,j >> GCA->log2s*((GCA->params->k-1)/2));
			}
		}
	}
}

void testCompress(int argc, char** argv)
{
	unsigned int i,j,ii,N,k,t,rule,sum,max,max_i,presum,postsum,min,min_i,lastmin_i[P];
	unsigned char** flags;
	unsigned int count, avglen,maxlen,prev,len;
	unsigned int counts[256];
	GraphCellularAutomaton **ECAs;
	chunk *ic;
	chunk **sols;
	chunk *preImages;
	unsigned int n;
	N = (unsigned int)atoi(argv[1]);
	k = (unsigned int)atoi(argv[2]);
	rule = (unsigned int)atoi(argv[3]);
	
	ECAs = (GraphCellularAutomaton **)malloc(256*sizeof(GraphCellularAutomaton *));
	flags = (unsigned char**)malloc(256*sizeof(unsigned char *));
	
	for (i=0;i<256;i++){
		ECAs[i] = CreateECA(N,k,i);
	}
	ic = (chunk *)malloc(ECAs[0]->size*sizeof(chunk));
	for (i=0;i<ECAs[0]->size;i++)
	{
		ic[i] = rand();
	}
	min = ECAs[0]->params->N;
	min_i = 0;
	
	for (j=0;j<P;j++)
		lastmin_i[j] = 0;
	for (i=0;i<256;i++)
	{
		for (j=0;j<ECAs[0]->size;j++)
		{
			ECAs[i]->config[j] = ic[j];
		}
	}
	for (t=0;t<N/8;t++)
	{
	for (i=0;i<256;i++)
	{
		
		printf("rule: %d\n",i);
		/*for (j=0;j<ECAs[i]->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECAs[i],j,0));
		}
		printf("\n");
		preImages = CAGetPreImages(ECAs[i],&n,NULL);
		for (j=0;j<n;j++)
		{
			for (k=0;k<ECAs[i]->size;k++)
				ECAs[i]->config[k] = (preImages+j*ECAs[i]->size)[k];
			for (k=0;k<ECAs[i]->params->N;k++)
			{
				printf("%u",GetCellStatePacked(ECAs[i],k,0));
			}
			printf(" -> ");
			CANextStep(ECAs[i]);
			for (k=0;k<ECAs[i]->params->N;k++)
			{
				printf("%u",GetCellStatePacked(ECAs[i],k,0));
			}
			printf("\n");
		}*/
		GetLeastSetBits(ECAs[i],ECAs[i]->config);
		/*for (j=0;j<ECAs[i]->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECAs[i],j,0));
		}
		printf("\n");*/
		CANextStep(ECAs[i]);
		/*for (j=0;j<ECAs[i]->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECAs[i],j,0));
		}
		printf("\n");*/
		if (memcmp(ECAs[i]->config,ic,ECAs[i]->size) != 0)
		{
	//		printf("error...\n");
		}
		else{
			register int sum;
			sum = 0;
			for (j=0;j<ECAs[i]->params->N;j++)
			{
				sum += GetCellStatePacked(ECAs[i],j,1);
			}
			if (i == lastmin_i[0])
			{
				if (sum < min)
				{
					min = sum;
					min_i = i;
				}
			}
			else
			{
				if (sum <= min)
				{
					min = sum;
					min_i = i;
				}
			}
		}
		
		
		
	}
		printf("best %d : rule %d\n",min,min_i);
		
		for (j=0;j<ECAs[min_i]->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECAs[min_i],j,1));
		}
		printf(" -> ");
		for (j=0;j<ECAs[min_i]->params->N;j++)
		{
			printf("%u",GetCellStatePacked(ECAs[min_i],j,0));
		}
		printf("\n");
		
		for (j=0;j<ECAs[0]->size;j++)
			ic[j] = ECAs[min_i]->st_pattern[1][j];
		for (i=0;i<256;i++)
		{
			for (j=0;j<ECAs[0]->size;j++)
			{
				ECAs[i]->config[j] = ic[j];
			}
		}
		lastmin_i[0] = min_i;
	}
	return;
	for (t=0;t<N/8;t++)
	{
		printf("t=%d\n",t);
	
	//ECA = CreateECA(N,k,rule);
	/*for(i=0;i<ECAs[110]->params->N;i++)
	{
		printf(" %u",GetCellStatePacked(ECAs[110],i,0));
	}
	printf("\n");*/
	//ECA->config[0] = 0x18f6;
		for (i=0;i<256;i++)
		{
			for (j=0;j<ECAs[0]->size;j++)
			{
				ECAs[i]->config[j] = ic[j];
			}
		}
	
		presum = 0;
		for (j=0;j<ECAs[0]->params->N;j++)
		{
			presum += GetCellStatePacked(ECAs[0],j,0);
		}
		min = presum;
		max = 0;
		count = 0;
		prev = GetCellStatePacked(ECAs[0],0,0);
		len = 0;
		maxlen = 0;
		for (j=1;j<ECAs[0]->params->N;j++)
		{
			if (prev != GetCellStatePacked(ECAs[0],j,0)) {
			count++;
				if (len > maxlen)
				{
					maxlen = len;
				}
				len = 0;
				prev = GetCellStatePacked(ECAs[0],j,0);
			}else{
				len++;
			}
		}
		min = presum;
		printf("Orignal set bits = %d\n",min);
		for (i=0;i<ECAs[0]->params->N;i++)
		{
			printf("%d",GetCellStatePacked(ECAs[0],i,0));
		}
		printf("\n");
		min += 5;
		for (i=0;i<256;i++)
		{
			unsigned char b;
			b = 0;
			for (j=0;j<P;j++)
				b |= (i == lastmin_i[j]);
			if(!b){
				//min_i = 110; 
				//printf("rule: %u",ECAs[i]->params->rule);
				/*for(j=0;j<ECAs[i]->params->N;j++)
			{
				printf(" %u",GetCellStatePacked(ECAs[i],j,0));
			}
			printf("\n");
			*/
			flags[i] = GetFlags(ECAs[i]);
			sum = 0;
			for (j=0;j<ECAs[i]->LUT_size*N;j++)
			{
				sum +=flags[i][j];
			}
		
			if (sum > 0)
			{	
			
			GetLeastBits(ECAs[i],flags[i]);
			postsum = 0;
			for (j=0;j<ECAs[i]->params->N;j++)
			{
				postsum += GetCellStatePacked(ECAs[i],j,0);
			}
			
			count = 0;
			len = 0;
			maxlen = 0;
			prev = GetCellStatePacked(ECAs[i],0,0);
			for (j=1;j<ECAs[0]->params->N;j++)
			{
				if (prev != GetCellStatePacked(ECAs[i],j,0)) {
					count++;
					if (len > maxlen)
					{
						maxlen = len;
					}
					len = 0;
					prev = GetCellStatePacked(ECAs[i],j,0);
				}else{
					len++;
				}
			}
			
			
			/*if (maxlen > max)
			{
				max = maxlen;
				max_i = i;
			}*/
			/*for(j=0;j<ECAs[i]->params->N;j++)
			{
				printf(" %u",GetCellStatePacked(ECAs[i],j,0));
			}
			printf("\n");
			*/
			/*test pre-image is correct: just testing for bugs*/
			CANextStep(ECAs[i]);
			
			/*for (j=0;j<ECAs[i]->LUT_size;j++)
			{
				for (ii=0;ii<ECAs[i]->params->N;ii++)
				{
					printf (" %d",flags[i][ii*ECAs[i]->LUT_size +j]);
				}
				printf("\n");
			}*/
			if (memcmp(ECAs[i]->config,ic,ECAs[0]->size) != 0)
			{
				//printf("error...\n");
				//exit(1);
			}
			else{
				counts[i] = count;
			if (postsum < min )
			{
				min = postsum;
				min_i = i;
				
			}
			}
		/*	for(j=0;j<ECAs[i]->params->N;j++)
			{
				printf(" %u",GetCellStatePacked(ECAs[i],j,0));
			}
			printf("\n");*/
		
		}
		else
		{
			counts[i] = RAND_MAX;
		}
		}
		/*else
		{
			printf("GOE\n");
		}*/
	}
	/*for (i=0;i<256;i++)
	{
		unsigned char b;
		b = 0;
		for (j=0;j<P;j++)
			b |= (i == lastmin_i[j]);
		if ((counts[i] <= min) && !b)
		{
			printf("get here? %d\n",min_i);
			min_i = i;
			break;
		}
	}*/
	
	/*count = 0;
			len = 0;
	
			prev = GetCellStatePacked(ECAs[i],0,0);
			for (j=1;j<ECAs[0]->params->N;j++)
			{
				if (prev != GetCellStatePacked(ECAs[max_i],j,1)) {
					if (len > 8)
					{
						count++;
					}
					len = 0;
					prev = GetCellStatePacked(ECAs[max_i],j,1);
				}else{
					len++;
				}
			}
	*/
	printf("Min set bits was %d using rule %d,%d\n",min,min_i,count);
	for (i=0;i<ECAs[0]->params->N;i++)
	{
		printf("%d",GetCellStatePacked(ECAs[min_i],i,0));
		//printf("%d",GetCellStatePacked(ECAs[max_i],i,0));
	
	}
	printf("\n");
	for (i=0;i<ECAs[0]->params->N;i++)
	{
		printf("%d",GetCellStatePacked(ECAs[min_i],i,1));
		//printf("%d",GetCellStatePacked(ECAs[max_i],i,1));
	}
	printf("\n");
	for (i=0;i<ECAs[0]->size;i++)
	{
		ic[i] = ECAs[min_i]->st_pattern[1][i];
//		ic[i] = ECAs[max_i]->st_pattern[1][i];

	}
	for (i=1;i<P;i++)
	{
		lastmin_i[i-1] = lastmin_i[i];
	}
		lastmin_i[P-1] = min_i;
	//flags = GetFlags(ECA);
	
	/*for (j=0;j<ECAs[110]->LUT_size;j++)
	{
		for (i=0;i<ECAs[110]->params->N;i++)
		{
			printf(" %d",flags[110][i*(ECAs[110]->LUT_size)+j]);	
		}
		printf(" %d \n",j);
	}
	printf("\n");*/
	}
	//GetLeastBits(ECAs[110],flags[110]);
	
	/*for (j=0;j<ECAs[110]->LUT_size;j++)
	{
		for (i=0;i<ECAs[110]->params->N;i++)
		{
			printf(" %d",flags[110][i*(ECAs[110]->LUT_size)+j]);	
		}
		printf(" %d \n",j);
	}*/
	
	
}

int testAttrCycles(int argc, char ** argv)
{
	GraphCellularAutomaton *ECA;
	chunk ic;
	chunk i;
	unsigned char* bin[8] = {"000","100","010","110","001","101","011","111"};
	unsigned int j,k;
	unsigned char* flags;
	ECA = CreateECA(30,3,30);
	printf("%d %d\n",ECA->params->N,ECA->LUT_size);
	for (i=0;i<16;i++)
	{
		ic = i;
		SetCAIC(ECA,&ic,EXPLICIT_IC_TYPE);
		flags = GetFlags(ECA);
		/*print the cell state*/
		printf("IC: %d\n",i);
		for (j=0;j<ECA->params->N;j++)
		{
			printf("%hhu",GetCellStatePacked(ECA,j,0));
		}
		printf("\n");
		printf("flags array\n");
		for(k=0;k<ECA->LUT_size;k++)
		{
			printf("%s : ",bin[k]);
			for (j=0;j<ECA->params->N;j++)
			{
				printf("%hhu ",flags[j*(ECA->LUT_size)+k]);
			}
			printf("\n");
		}
	}
}

int main(int argc, char** argv)
{
	/*testbitaccess(argc,argv);*/
	/*testECA(argc,argv);*/
	/*testRevAlgorithm(argc,argv);*/
	/*testCompress(argc,argv);*/
	/*testSE(argc,argv);*/
	testAttrCycles(argc,argv);
}
