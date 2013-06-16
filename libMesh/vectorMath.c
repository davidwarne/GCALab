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
/* File: vectorMath.c
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 10/03/2012
 * Last Modified: 19/03/2012
 *
 * Version History:
 *       v 0.01 (10/03/2012) - Initial Version
 *       v 0.02 (10/03/2012) - Implemented Add_f,Diff_f,MidPoint_f,Mean_f,Norm_f
 *                             NormDiff_f,Dot_f,Cross_f,MagCross_f,ScalarTriple_f,
 *                             Add_lf,Diff_lf,MidPoint_lf,Mean_lf,Norm_lf
 *                             NormDiff_lf,Dot_lf,Cross_lf,MagCross_lf,ScalarTriple_lf
 *       v 0.022 (19/03/2012) - added Normal_f() and Normal_lf()
 *
 * Descritpion: library of common vector math functions
 *
 * TODO:
 *    Nothin'
 *
 * Known Issues:
 *     There are currently no known issues
 * =============================================================================
 */

#include "vectorMath.h"

/* Add_f(): computes the vector addition of vectors x,y (32-bit floats) in R^n
 *
 * Parameters:
 *     x,y - vectors to compute x + y
 *     n - dimension of x and y
 * Returns:
 *     z - the sum of x and y
 */
float * Add_f(float *x, float *y, int n)
{
	float * z;
	int i;
	
	if (!(z = (float *)malloc(n*sizeof(float))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

/* Diff_f(): computes the vector difference of vectors x,y (32-bit floats) in R^n
 *
 * Parameters:
 *     x,y - vectors to compute x - y
 *     n - dimension of x and y
 * Returns:
 *     z - the difference of x and y
 */
float * Diff_f(float *x, float *y, int n)
{
	float * z;
	int i;
	
	if (!(z = (float *)malloc(n*sizeof(float))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] - y[i];
	}
	return z;
}

/* MidPoint_f(): calculates the midpoint of two single vectors (32-bit float) 
 *               of dimension n
 * Parameters:
 *     x,y - two vectors to compute midpoint of
 *     n - dimension of x and y
 * Returns:
 *     m - pointer to memory containin the midpoint
 * Pre-condition:
 *	   n must be valid for both x and y
 */
float *Midpoint_f(float *x, float *y, int n)
{
	float *m;
	int i;
	
	if (!(m = (float *)malloc(n*sizeof(float))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		m[i] = (x[i] + y[i])*0.5;
	}
	return m;
}

/* Mean_f(): computes the mean of num vectors in R^n
 *
 * Parameters:
 *     num - number of vectors to average
 *     vList - list of vectors
 *     n - dimension of the vectors
 * Returns:
 *     the mean vector
 * Pre-condition:
 *     num must be equal to the number of vectors in the variab;e argument list
 */
float *Mean_f(int num,float * vList, int n)
{
	float *m,*vec;
	int i,j;
	if (!(m = (float *)malloc(n*sizeof(float))))
	{
		return NULL;
	}
	
	for (j=0;j<n;j++)
	{
		m[j] = vList[j];
	}
	
	for (i=1;i<num;i++)
	{
		vec = vList + i*n;
		for (j=0;j<n;j++)
		{
			m[j] += vec[j];
		}
	}
	
	for (j=0;j<n;j++)
	{
		m[j] /= (float)n;
	}
	
	return m;
}

/* Norm_f(): calculates the norm of a vector (32-bit float) of dimension n
 * 
 * Parameters:
 *     x - vector to compute norm of
 *     n - dimension of x
 *     norm - the type of norm to take support 1-norm, 2-norm or infinity-norm
 * Returns:
 *     n - the norm of x
 * Pre-condition:
 *	   nrm must be valid for x
 * Note: 
 *     x and y are interpreted as column vectors for the purpose of the 1-norm
 *     or the infinity-norm
 */
float Norm_f(float *x,int n,unsigned char norm)
{
	float nrm;
	int i;
	
	nrm = 0.0;
	switch(norm)
	{
		/*the sum of the elements*/
		case ONE_NORM:
			for (i=0;i<n;i++)
			{
				nrm += fabs(x[i]);
			}
			break;
		/*max element*/
		case INF_NORM:
			for (i=0;i<n;i++)
			{
				if (nrm < fabs(x[i]))
				{
					nrm = fabs(x[i]);
				}
			}
			break;
		/* euclidean norm - sqrt of sum of squares*/
		default:
		case TWO_NORM:
			for (i=0;i<n;i++)
			{
				nrm += x[i]*x[i];
			}
			nrm = sqrt(nrm);
			break;
	}
	return nrm;
}

/* NormDiff_f(): calculates the norm of the differenc of two vectors (32-bit float)
 *               of dimension n
 * Parameters:
 *     x,y - two vectors to compute norm diff
 *     n - dimension of x and y
 *     norm - the type of norm to take support 1-norm, 2-norm or infinity-norm
 * Returns:
 *     nd - the norm of the difference
 * Pre-condition:
 *	   n must be valid for both x and y
 * Note: 
 *     x and y are interpreted as column vectors for the purpose of the 1-norm
 *     or the infinity-norm
 */
float NormDiff_f(float *x, float *y,int n,unsigned char norm)
{
	float nrm, *z;
	int i;
	
	nrm = 0.0;
	z = Diff_f(x,y,n);
	nrm = Norm_f(z,n,norm);
	free(z);
	return nrm;
}

/* Dot_f(): vector inner product in R^n
 *
 * Parameters:
 *     x,y - two vectors in R^n
 *     n - dimension of x and y
 * Return:
 *     vector inner product (a scalar)
 * Pre-condition:
 *	   n must be valid for both x and y
 */
float Dot_f(float *x, float *y, int n)
{
	float d;
	int i;
	d = 0;
	for (i=0;i<n;i++)
	{
		d += x[i]*y[i];
	}
	return d;
}

/* Cross_f(): vector cross product in R^3
 *
 * Parameters:
 *     x,y - two vectors in R^3
 * Return:
 *     vector cross product (a vector)
 * Pre-condition:
 *	   both x and y must be in R^3
 */
float *Cross_f(float *x, float *y)
{
	float *z;
	if (!(z = (float *)malloc(3*sizeof(float))))
	{
		return NULL;
	}
	
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
	
	return z;
}

/* MagCross_f(): magnitude vector cross product in R^3
 *
 * Parameters:
 *     x,y - two vectors in R^3
 * Return:
 *     magnitude of vector cross product (a scalar)
 * Pre-condition:
 *	   both x and y must be in R^3
 */
float MagCross_f(float *x, float *y)
{
	float z0,z1,z2;
	z0 = x[1]*y[2] - x[2]*y[1];
	z1 = x[2]*y[0] - x[0]*y[2];
	z2 = x[0]*y[1] - x[1]*y[0];

	return sqrt(z0*z0+z1*z1+z2*z2);
}

/* Normal_f(): computes the normal to the plane defined by the three given poitns
 *
 * Parameters:
 *     x,y,z - vectors in R^3
 */
float *Normal_f(float *x, float *y,float *z)
{
	float *v0;
	float *v1;
	float *normal;
	float mag;
	v0 = Diff_f(y,x,3);
	v1 = Diff_f(z,x,3);
	
	normal = Cross_f(v0,v1);
	mag = Norm_f(normal,3,TWO_NORM);
	normal[0] /= mag;
	normal[1] /= mag;
	normal[2] /= mag;
	free(v0);
	free(v1);
	return normal;
}

/* ScalarTriple_f(): scalar triple product in R^3 
 *
 * Parameters:
 *     x,y,z - three vectors in R^3
 * Return:
 *     dot(x,cross(y,z))
 * Pre-condition:
 *	   x,y and z must be in R^3
 */
float ScalarTriple_f(float *x, float *y,float *z)
{
	float *tmp,stp;
	tmp = Cross_f(y,z);
	stp = Dot_f(x,tmp,3);
	free(tmp);
	return stp;
}

/* Add_lf(): computes the vector addition of vectors x,y (64-bit floats) in R^n
 *
 * Parameters:
 *     x,y - vectors to compute x + y
 *     n - dimension of x and y
 * Returns:
 *     z - the sum of x and y
 */
double * Add_lf(double *x, double *y, int n)
{
	double * z;
	int i;
	
	if (!(z = (double *)malloc(n*sizeof(double))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

/* Diff_lf(): computes the vector difference of vectors x,y (64-bit floats) in R^n
 *
 * Parameters:
 *     x,y - vectors to compute x - y
 *     n - dimension of x and y
 * Returns:
 *     z - the difference of x and y
 */
double * Diff_lf(double *x, double *y, int n)
{
	double * z;
	int i;
	
	if (!(z = (double *)malloc(n*sizeof(double))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] - y[i];
	}
	return z;
}

/* MidPoint_lf(): calculates the midpoint of two single vectors (64-bit float) 
 *               of dimension n
 * Parameters:
 *     x,y - two vectors to compute midpoint of
 *     n - dimension of x and y
 * Returns:
 *     m - pointer to memory containin the midpoint
 * Pre-condition:
 *	   n must be valid for both x and y
 */
double *Midpoint_lf(double *x, double *y, int n)
{
	double *m;
	int i;
	
	if (!(m = (double *)malloc(n*sizeof(double))))
	{
		return NULL;
	}
	
	for (i=0;i<n;i++)
	{
		m[i] = (x[i] + y[i])*0.5;
	}
	return m;
}

/* Mean_lf(): computes the mean of num vectors in R^n
 *
 * Parameters:
 *     num - number of vectors to average
 *     vList - list of vectors
 *     n - dimension of the vectors
 * Returns:
 *     the mean vector
 * Pre-condition:
 *     num must be equal to the number of vectors in the variab;e argument list
 */
double *Mean_lf(int num,double * vList, int n)
{
	double *m,*vec;
	int i,j;
	if (!(m = (double *)malloc(n*sizeof(double))))
	{
		return NULL;
	}
	
	for (j=0;j<n;j++)
	{
		m[j] = vList[j];
	}
	
	for (i=1;i<num;i++)
	{
		vec = vList + i*n;
		for (j=0;j<n;j++)
		{
			m[j] += vec[j];
		}
	}
	
	for (j=0;j<n;j++)
	{
		m[j] /= (double)n;
	}
	
	return m;
}

/* Norm_lf(): calculates the norm of a vector (64-bit float) of dimension n
 * 
 * Parameters:
 *     x - vector to compute norm of
 *     n - dimension of x
 *     norm - the type of norm to take support 1-norm, 2-norm or infinity-norm
 * Returns:
 *     n - the norm of x
 * Pre-condition:
 *	   nrm must be valid for x
 * Note: 
 *     x and y are interpreted as column vectors for the purpose of the 1-norm
 *     or the infinity-norm
 */
double Norm_lf(double *x,int n,unsigned char norm)
{
	double nrm;
	int i;
	
	nrm = 0.0;
	switch(norm)
	{
		/*the sum of the elements*/
		case ONE_NORM:
			for (i=0;i<n;i++)
			{
				nrm += x[i];
			}
			break;
		/*max element*/
		case INF_NORM:
			for (i=0;i<n;i++)
			{
				if (nrm < x[i])
				{
					nrm = x[i];
				}
			}
			break;
		/* euclidean norm - sqrt of sum of squares*/
		default:
		case TWO_NORM:
			for (i=0;i<n;i++)
			{
				nrm += x[i]*x[i];
			}
			nrm = sqrt(nrm);
			break;
	}
	return nrm;
}

/* NormDiff_lf(): calculates the norm of the differenc of two vectors (64-bit float)
 *               of dimension n
 * Parameters:
 *     x,y - two vectors to compute norm diff
 *     n - dimension of x and y
 *     norm - the type of norm to take support 1-norm, 2-norm or infinity-norm
 * Returns:
 *     nd - the norm of the difference
 * Pre-condition:
 *	   n must be valid for both x and y
 * Note: 
 *     x and y are interpreted as column vectors for the purpose of the 1-norm
 *     or the infinity-norm
 */
double NormDiff_lf(double *x, double *y,int n,unsigned char norm)
{
	double nrm, *z;
	int i;
	
	nrm = 0.0;
	z = Diff_lf(x,y,n);
	Norm_lf(z,n,norm);
	free(z);
	return nrm;
}

/* Dot_lf(): vector inner product in R^n
 *
 * Parameters:
 *     x,y - two vectors in R^n
 *     n - dimension of x and y
 * Return:
 *     vector inner product (a scalar)
 * Pre-condition:
 *	   n must be valid for both x and y
 */
double Dot_lf(double *x, double *y, int n)
{
	double d;
	int i;
	d = 0;
	for (i=0;i<n;i++)
	{
		d += x[i]*y[i];
	}
	return d;
}

/* Cross_lf(): vector cross product in R^3
 *
 * Parameters:
 *     x,y - two vectors in R^3
 * Return:
 *     vector cross product (a vector)
 * Pre-condition:
 *	   both x and y must be in R^3
 */
double *Cross_lf(double *x, double *y)
{
	double *z;
	if (!(z = (double *)malloc(3*sizeof(double))))
	{
		return NULL;
	}
	
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
	
	return z;
}

/* MagCross_lf(): magnitude vector cross product in R^3
 *
 * Parameters:
 *     x,y - two vectors in R^3
 * Return:
 *     magnitude of vector cross product (a scalar)
 * Pre-condition:
 *	   both x and y must be in R^3
 */
double MagCross_lf(double *x, double *y)
{
	double z0,z1,z2;
	z0 = x[1]*y[2] - x[2]*y[1];
	z1 = x[2]*y[0] - x[0]*y[2];
	z2 = x[0]*y[1] - x[1]*y[0];

	return sqrt(z0*z0+z1*z1+z2*z2);
}

/* Normal_lf(): computes the normal to the plane defined by the three given poitns
 *
 * Parameters:
 *     x,y,z - vectors in R^3
 */
double *Normal_lf(double *x, double *y,double *z)
{
	double *v0;
	double *v1;
	double *normal;
	v0 = Diff_lf(y,x,3);
	v1 = Diff_lf(z,x,3);
	
	normal = Cross_lf(v0,v1); 
	free(v0);
	free(v1);
	return normal;
}

/* ScalarTriple_f(): scalar triple product in R^3 
 *
 * Parameters:
 *     x,y,z - three vectors in R^3
 * Return:
 *     dot(x,cross(y,z))
 * Pre-condition:
 *	   x,y and z must be in R^3
 */
double ScalarTriple_lf(double *x, double *y,double *z)
{
	double *tmp,stp;
	tmp = Cross_lf(y,z);
	stp = Dot_lf(x,tmp,3);
	free(tmp);
	return stp;
}


