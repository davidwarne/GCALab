/* File: GCALab_fio.c
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 19/04/2012
 * Last Modified: 19/04/2012
 *
 * Description: File IO function implementations for GCALab formats
 *
 *==============================================================================
 */
 
#include "GCALab_fio.h"

char *GCALab_fio_format[14] = {"%f","%lf","%hhu","%hu","%u","%llu","%hhd","%hd","%d","%lld","%hhx","%hx","%x","%llx"};

/* GCALAB_saveCA(): Writes CA data to a file
 *
 * Parameters:
 *     filename - the file to write to
 *     GCA - Graph Cellular Automaton to write to file
 *
 * Note:
 *   This function will NOT write any spatio-temporal information, only sufficient
 *   information to define the rule.
 */ 
char GCALab_fio_saveCA(char *filename,GraphCellularAutomaton *GCA,mesh *m)
{
	FILE* fp;
	unsigned int i,j;
	char lutfile[255],meshfile[255],graphfile[255];
	if (!(fp = fopen(filename,"w")))
	{
		return WRITE_FAILED;
	}
	
	sprintf(lutfile,"%s.lut",filename);
	sprintf(meshfile,"%s.off",filename);
	sprintf(graphfile,"%s.top",filename);

	fprintf(fp,"#\n# Created by GCALab.\n#\n");
	fprintf(fp,"#\n#CA Parameters\n#\n");
	fprintf(fp,"N : %d\n",GCA->params->N);
	fprintf(fp,"s : %u\n",GCA->params->s);
	fprintf(fp,"k : %u\n",GCA->params->k);
	fprintf(fp,"type : %u\n",GCA->params->rule_type);
	fprintf(fp,"rule : %u\n",GCA->params->rule);
	fprintf(fp,"WSIZE : %d\n",GCA->params->WSIZE);
	fprintf(fp,"LUT : %s\n",lutfile);
	fprintf(fp,"topology : %s\n",graphfile);
	if (m != NULL)
	{
		fprintf(fp,"mesh : %s\n",meshfile);
	}
	fclose(fp);

    	
	if (!(fp = fopen(lutfile,"w")))
	{
		return WRITE_FAILED;
	}

	for (i=0;i<GCA->LUT_size;i++)
	{
		fprintf(fp,"%u,%u\n",i,GCA->ruleLUT[i]);
	}
	fclose(fp);
	

	if (!(fp = fopen(graphfile,"w")))
	{
		return WRITE_FAILED;
	}
	for (i=0;i<GCA->params->N;i++)
	{
		fprintf(fp,"%d",i);
		for (j=0;j<GCA->params->k-1;j++)
		{
			fprintf(fp,",%d",GCA->params->graph[i*(GCA->params->k-1) +j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	if (m != NULL)
	{
		SaveMesh(meshfile,m,OFF_FORMAT);
	}
	return WRITE_SUCCESS;
};

char GCALab_fio_saveData(char* filename,char * name, void * data, int N,unsigned char type)
{
	FILE* fp;
	int i;

	if (!(fp = fopen(filename,"w")))
	{
		return WRITE_FAILED;
	}
	
	switch(type)
	{
		case FLOAT32:
		{
			float *data_fp32;
			data_fp32 = (float*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %f",data_fp32[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case FLOAT64:
		{
			double *data_fp64;
			data_fp64 = (double*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %lf",data_fp64[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case UINT8:
		{
			unsigned char *data_uint8;
			data_uint8 = (unsigned char*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %hhu",data_uint8[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case UINT16:
		{
			unsigned short *data_uint16;
			data_uint16 = (unsigned short*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %hu",data_uint16[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case UINT32:
		{
			unsigned int *data_uint32;
			data_uint32 = (unsigned int*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %u",data_uint32[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case UINT64:
		{
			unsigned long long *data_uint64;
			data_uint64 = (unsigned long long*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %llu",data_uint64[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case SINT8:
		{
			signed char *data_uint8;
			data_uint8 = (signed char*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %hhd",data_uint8[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case SINT16:
		{
			signed short *data_uint16;
			data_uint16 = (signed short*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %hd",data_uint16[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case SINT32:
		{
			signed int *data_uint32;
			data_uint32 = (signed int*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %d",data_uint32[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case SINT64:
		{
			signed long long *data_uint64;
			data_uint64 = (signed long long*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %lld",data_uint64[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case HEX8:
		{
			unsigned char *data_uint8;
			data_uint8 = (unsigned char*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %02hhx",data_uint8[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case HEX16:
		{
			unsigned short *data_uint16;
			data_uint16 = (unsigned short*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %04hx",data_uint16[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case HEX32:
		{
			unsigned int *data_uint32;
			data_uint32 = (unsigned int*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %08x",data_uint32[i]);
			}
			fprintf(fp,"\n");
		}
			break;
		case HEX64:
		{
			unsigned long long *data_uint64;
			data_uint64 = (unsigned long long*)data;
			for (i=0;i<N;i++)
			{
				fprintf(fp," %016llx",data_uint64[i]);
			}
			fprintf(fp,"\n");
		}
			break;
	}
	fclose(fp);
	return WRITE_SUCCESS;
}
