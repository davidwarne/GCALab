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
 *       v 0.03 (19/06/2013) - Finally go round to fixinf the  ssue with many internal
 *                             mallocs in the functions. Now optional memory may be supplied by the calling funtion.
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

/** @brief Computes the vector addition of vectors x,y (32-bit floats) in R^n
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param z_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the element wise sum of \a x and \a y.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
float * Add_f(float *x, float *y,float* z_out,int n)
{
	float * z;
	int i;
	
    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (float *)malloc(n*sizeof(float))))
	    {
		    return NULL;
	    }
    }
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

/** @brief Computes the vector difference of vectors x,y (32-bit floats) in R^n
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param z_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the element wise differnce of \a x and \a y.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
float * Diff_f(float *x, float *y, float* z_out,int n)
{
	float * z;
	int i;
	
    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (float *)malloc(n*sizeof(float))))
	    {
		    return NULL;
	    }
	}

	for (i=0;i<n;i++)
	{
		z[i] = x[i] - y[i];
	}
	return z;
}

/** @brief Calculates the midpoint of two single vectors (32-bit float) 
 *  of dimension n.
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param m_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the midpoint.
 * @note If \a m_out != NULL then the returned point is equal to \a m_out. Otherwise
 * new memory is allocted.
 */
float *Midpoint_f(float *x, float *y, float *m_out, int n)
{
	float *m;
	int i;
	
    if (m_out != NULL)
    {
        m = m_out;
    }
    else
    {
	    if (!(m = (float *)malloc(n*sizeof(float))))
	    {
		    return NULL;
	    }
    }
	
	for (i=0;i<n;i++)
	{
		m[i] = (x[i] + y[i])*0.5;
	}
	return m;
}

/** @brief Computes the mean of num vectors in R^n.
 *
 * @param num The number of vectors to average.
 * @param vList A list of vectors.
 * @param m_out User supplied pointer for results.
 * @param n The dimension of the vectors.
 * 
 * @returns A pointer to the mean vector.
 * @note If \a m_out != NULL then the returned point is equal to \a m_out. Otherwise
 * new memory is allocted.
 */
float *Mean_f(int num,float * vList, float* m_out, int n)
{
	float *m,*vec;
	int i,j;
    if (m_out != NULL)
    {
        m = m_out;
    }
    else
    {
	    if (!(m = (float *)malloc(n*sizeof(float))))
	    {
		    return NULL;
	    }
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

/** @brief Calculates the norm of a vector (32-bit float) of dimension n.
 * 
 * @param x Vector to compute norm of.
 * @param n Dimension of x.
 * @param norm The type of norm to take support 1-norm, 2-norm or infinity-norm.
 *
 * @returns The norm of x.
 *
 * @note \a x and \a y are interpreted as column vectors for the purpose of the 1-norm
 * or the infinity-norm
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

/** @brief Calculates the norm of the difference of two vectors (32-bit float)
 *  of dimension n
 *
 * @param x First vector operand. 
 * @param y Second vector operand.
 * @param tmp User supplied memory to avoid repeated mallocs.
 * @param n Dimension of \a x and \a y.
 * @param norm The type of norm to take support 1-norm, 2-norm or infinity-norm.
 *
 * @returns The norm of the difference.
 * @note \a x and \a y are interpreted as column vectors for the purpose of the 1-norm
 * or the infinity-norm.
 * @note If \a tmp == NULL the function will still operate correctly, but will perform
 * internal mallocs.
 */
float NormDiff_f(float *x, float *y,float* tmp,int n,unsigned char norm)
{
	float nrm, *z;
	int i;
	
    z = tmp;
	nrm = 0.0;
	z = Diff_f(x,y,z,n);
	nrm = Norm_f(z,n,norm);
	return nrm;
}

/** @brief Vector inner product in R^n
 *
 * @param x First vector operand.
 * @param y Second vectot operand.
 * @param n Dimension of \a x and \a y. 
 * 
 * @returns The vector inner product (a scalar).
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

/** @brief The vector cross product in R^3
 *
 * @param x The first vector operand in R^3.
 * @param y The second vector operand in R^3.
 * @param z_out User supplied memory for result.
 *
 * @returns vector cross product (a vector).
 * @remark Both x and y must be in R^3.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
float *Cross_f(float *x, float *y,float *z_out)
{
	float *z;

    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (float *)malloc(3*sizeof(float))))
	    {
		    return NULL;
	    }
    }
	
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
	
	return z;
}

/** @brief Computes magnitude vector cross product in R^3.
 *
 * @param x The first vector operand in R^3.
 * @param y The second vector operand in R^3.
 * 
 * @returns The magnitude of vector cross product (a scalar).
 * @remark Both x and y must be in R^3.
 * @remark The magnitude is equal to the area of the parallelogram subtended by \a x and \a y.
 */
float MagCross_f(float *x, float *y)
{
	float z0,z1,z2;
	z0 = x[1]*y[2] - x[2]*y[1];
	z1 = x[2]*y[0] - x[0]*y[2];
	z2 = x[0]*y[1] - x[1]*y[0];

	return sqrt(z0*z0+z1*z1+z2*z2);
}

/** @brief computes the unit normal to the plane defined by the three given points.
 *
 * @param x The first point operand.
 * @param y The second point operand.
 * @param z The third point operand.
 * @param n_out User suplied memory for the unit normal output.
 *
 * @returns A pointer to the unit normal to the plain defined by 
 * \a x \a y and \a z.
 *
 * @note If \a n_out != NULL then the returned point is equal to \a n_out. Otherwise new memory is allocated.
 */
float *Normal_f(float *x, float *y,float *z, float* n_out)
{
	float v0[3];
	float v1[3];
	float *normal;
	float inv_mag;

	normal = n_out;

    /*Get two linearly independent vectors*/
	Diff_f(y,x,v0,3);
	Diff_f(z,x,v1,3);
	/*compute the unit normal*/
	normal = Cross_f(v0,v1,normal);
	inv_mag = 1.0/Norm_f(normal,3,TWO_NORM);
	normal[0] *= inv_mag;
	normal[1] *= inv_mag;
	normal[2] *= inv_mag;

	return normal;
}

/** @brief The scalar triple product in R^3. 
 *
 * @param x The first point operand.
 * @param y The second point operand.
 * @param z The third point operand.
 *
 * @returns The scalar triple product defined by dot(x,cross(y,z))
 * @remark The scalar triple product gives the volume of the 
 * parallelpiped subtended by the vectors \a x \a y and \a z.
 */
float ScalarTriple_f(float *x, float *y,float *z)
{
	float tmp[3];
	float stp;
	Cross_f(y,z,tmp);
	stp = Dot_f(x,tmp,3);
	return stp;
}

/** @brief Computes the vector addition of vectors x,y (64-bit floats) in R^n
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param z_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the element wise sum of \a x and \a y.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
double * Add_lf(double *x, double *y,double* z_out,int n)
{
	double * z;
	int i;
	
    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (double *)malloc(n*sizeof(double))))
	    {
		    return NULL;
	    }
    }
	
	for (i=0;i<n;i++)
	{
		z[i] = x[i] + y[i];
	}
	return z;
}

/** @brief Computes the vector difference of vectors x,y (64-bit floats) in R^n
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param z_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the element wise differnce of \a x and \a y.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
double * Diff_lf(double *x, double *y, double* z_out,int n)
{
	double * z;
	int i;
	
    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (double *)malloc(n*sizeof(double))))
	    {
		    return NULL;
	    }
	}

	for (i=0;i<n;i++)
	{
		z[i] = x[i] - y[i];
	}
	return z;
}

/** @brief Calculates the midpoint of two single vectors (64-bit float) 
 *  of dimension n.
 *
 * @param x Fist vector operand.
 * @param y Second vector operand.
 * @param m_out user supplied pointer for results
 * @param n Dimension of operands.
 *
 * @returns pointer to vector which is the midpoint.
 * @note If \a m_out != NULL then the returned point is equal to \a m_out. Otherwise
 * new memory is allocted.
 */
double *Midpoint_lf(double *x, double *y, double *m_out, int n)
{
	double *m;
	int i;
	
    if (m_out != NULL)
    {
        m = m_out;
    }
    else
    {
	    if (!(m = (double *)malloc(n*sizeof(double))))
	    {
		    return NULL;
	    }
    }
	
	for (i=0;i<n;i++)
	{
		m[i] = (x[i] + y[i])*0.5;
	}
	return m;
}

/** @brief Computes the mean of num vectors in R^n.
 *
 * @param num The number of vectors to average.
 * @param vList A list of vectors.
 * @param m_out User supplied pointer for results.
 * @param n The dimension of the vectors.
 * 
 * @returns A pointer to the mean vector.
 * @note If \a m_out != NULL then the returned point is equal to \a m_out. Otherwise
 * new memory is allocted.
 */
double *Mean_lf(int num,double * vList, double* m_out, int n)
{
	double *m,*vec;
	int i,j;
    if (m_out != NULL)
    {
        m = m_out;
    }
    else
    {
	    if (!(m = (double *)malloc(n*sizeof(double))))
	    {
		    return NULL;
	    }
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

/** @brief Calculates the norm of a vector (64-bit float) of dimension n.
 * 
 * @param x Vector to compute norm of.
 * @param n Dimension of x.
 * @param norm The type of norm to take support 1-norm, 2-norm or infinity-norm.
 *
 * @returns The norm of x.
 *
 * @note \a x and \a y are interpreted as column vectors for the purpose of the 1-norm
 * or the infinity-norm
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

/** @brief Calculates the norm of the difference of two vectors (64-bit float)
 *  of dimension n
 *
 * @param x First vector operand. 
 * @param y Second vector operand.
 * @param tmp User supplied memory to avoid repeated mallocs.
 * @param n Dimension of \a x and \a y.
 * @param norm The type of norm to take support 1-norm, 2-norm or infinity-norm.
 *
 * @returns The norm of the difference.
 * @note \a x and \a y are interpreted as column vectors for the purpose of the 1-norm
 * or the infinity-norm.
 * @note If \a tmp == NULL the function will still operate correctly, but will perform
 * internal mallocs.
 */
double NormDiff_lf(double *x, double *y,double* tmp,int n,unsigned char norm)
{
	double nrm, *z;
	int i;
	
    z = tmp;
	nrm = 0.0;
	z = Diff_lf(x,y,z,n);
	nrm = Norm_lf(z,n,norm);
	free(z);
	return nrm;
}

/** @brief Vector inner product in R^n
 *
 * @param x First vector operand.
 * @param y Second vectot operand.
 * @param n Dimension of \a x and \a y. 
 * 
 * @returns The vector inner product (a scalar).
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

/** @brief The vector cross product in R^3
 *
 * @param x The first vector operand in R^3.
 * @param y The second vector operand in R^3.
 * @param z_out User supplied memory for result.
 *
 * @returns vector cross product (a vector).
 * @remark Both x and y must be in R^3.
 * @note If \a z_out != NULL then the returned point is equal to \a z_out. Otherwise
 * new memory is allocted.
 */
double *Cross_lf(double *x, double *y,double *z_out)
{
	double *z;

    if (z_out != NULL)
    {
        z = z_out;
    }
    else
    {
	    if (!(z = (double *)malloc(3*sizeof(double))))
	    {
		    return NULL;
	    }
    }
	
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
	
	return z;
}

/** @brief Computes magnitude vector cross product in R^3.
 *
 * @param x The first vector operand in R^3.
 * @param y The second vector operand in R^3.
 * 
 * @returns The magnitude of vector cross product (a scalar).
 * @remark Both x and y must be in R^3.
 * @remark The magnitude is equal to the area of the parallelogram subtended by \a x and \a y.
 */
double MagCross_lf(double *x, double *y)
{
	double z0,z1,z2;
	z0 = x[1]*y[2] - x[2]*y[1];
	z1 = x[2]*y[0] - x[0]*y[2];
	z2 = x[0]*y[1] - x[1]*y[0];

	return sqrt(z0*z0+z1*z1+z2*z2);
}

/** @brief computes the unit normal to the plane defined by the three given points.
 *
 * @param x The first point operand.
 * @param y The second point operand.
 * @param z The third point operand.
 * @param n_out User suplied memory for the unit normal output.
 *
 * @returns A pointer to the unit normal to the plain defined by 
 * \a x \a y and \a z.
 *
 * @note If \a n_out != NULL then the returned point is equal to \a n_out. Otherwise new memory is allocated.
 */
double *Normal_lf(double *x, double *y,double *z, double* n_out)
{
	double v0[3];
	double v1[3];
	double *normal;
	double inv_mag;

	normal = n_out;

    /*Get two linearly independent vectors*/
	Diff_lf(y,x,v0,3);
	Diff_lf(z,x,v1,3);
	/*compute the unit normal*/
	normal = Cross_lf(v0,v1,normal);
	inv_mag = 1.0/Norm_lf(normal,3,TWO_NORM);
	normal[0] *= inv_mag;
	normal[1] *= inv_mag;
	normal[2] *= inv_mag;

	return normal;
}

/** @brief The scalar triple product in R^3. 
 *
 * @param x The first point operand.
 * @param y The second point operand.
 * @param z The third point operand.
 *
 * @returns The scalar triple product defined by dot(x,cross(y,z))
 * @remark The scalar triple product gives the volume of the 
 * parallelpiped subtended by the vectors \a x \a y and \a z.
 */
double ScalarTriple_lf(double *x, double *y,double *z)
{
	double tmp[3];
	double stp;
	Cross_lf(y,z,tmp);
	stp = Dot_lf(x,tmp,3);
	return stp;
}
