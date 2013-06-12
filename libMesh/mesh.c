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
/* File: mesh.h
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
 *       v 0.01 (10/03/2012) - Initial Version... my 6-month old son helped daddy
 *                             with the header file
 *       v 0.02 (10/03/2012) - Implemented CreateVertexList,InsertVertex,DeleteVertex
 *                             GetVertex_ptr,GetVertex_cp,CreateEdgeList,InsertEdge,
 *                             DeleteEdge,GetEdge_ptr,GetEdge_cp,EdgeMidPoint
 *                             EdgeLength,CreateFaceList,InsertFace,DeleteFace,
 *                             GetFace_ptr,GetFace_cp,CheckErr,PrintErrorMsg
 *       v 0.03 (13/03/2012) - Implemented FaceMdiPoint, FaceArea, CreateMesh,
 *                             CopyMesh,SrufaceArea,GeometricCentre
 *       v 0.035 (19/03/2012) - Implemented SaveMesh() and working on SubDivide()
 *       v 0.040 (05/08/2012) - I have neglected this code... time to get back into it
 *                             - Implemented Readers and writers
 *                             - performed some testing on SubDivide(), starting to 
 *                               feel that I've over complecated things
 *       v 0.045 (11/08/2012) - decided to remove the edge list, and all references to it
 *       v 0.050 (14/08/2012) - Re-Implemented SubDivide (TODO: 4), and implemented
 *                              RemoveDuplicateVertices().
 *       v 0.055 (28/09/2012) - Implemented CreateTopology functions, and CreateDual()
 *       v 0.060 (13/02/2013) - Added support for genus-2 topolog mesh creation.
 *
 *
 * Descritpion: Implementation of mesh construction and manipulation functions
 *
 * TODO: 
 *    1. finish implemeting face functions - done (v 0.03)
 *    2. finish implementing mesh functions
 *    3. Running into problems with Writers and Subdivide, I am starting to get annoyed
 *       with the edge list... its very awkward to maintain windings. I question if its
 *       needed, I may decide to have a face an vertex list only. - done (v 0.045)
 *    4. SubDivide() needs to be re-implemented as edge list no longer exists. - done (v 0.050)
 *
 * Known Issues:
 *     There are currently no known issues
 * =============================================================================
 */
 
#include "mesh.h"

/** @brief Allocates memory for a vertex list.
 *
 * @param size The number of vertices.
 * @param dim The dimension of vertices.
 *
 * @returns A pointer to allocated vertexList struct.
 * @retVal OUT_OF_MEMORY on failure.
 */
vertexList * CreateVertexList(int size,int dim)
{
	vertexList *vList;
	int n;
	
	if (!(vList = (vertexList *)malloc(sizeof(vertexList))))
	{
		return OUT_OF_MEMORY;
	}
	
	vList->numVerts = 0;
	vList->size = size;
	vList->dim = dim;
	
	n = size*dim;
	if (!(vList->verts = (float *)malloc(n*sizeof(float))))
	{
		return OUT_OF_MEMORY;
	}
	memset((void *)(vList->verts),0,n*sizeof(float));
	return vList;
}

/** @brief Adds a vertes to the given vertex list.
 *
 * @param vList A pointer to vertex list to insert into.
 * @param v The vertex to be inserted.
 * 
 * @retVal INSERT_SUCCESS Successful insertion.
 * @retVal INSERT_FAILED Indicates a failure error code.
 *
 * @note length of \a v == <em>vList->dim</em>.
 */
r_code InsertVertex(vertexList *vList,float *v)
{
	int i;
	if (vList->size <= vList->numVerts)
	{
		int newsize;
		float *newverts;
		newsize = vList->size + REALLOC_EXTRA;
		newverts = (float*)realloc(vList->verts,newsize*(vList->dim)*sizeof(float));
		if (!newverts)
		{
			return INSERT_FAILED;
		}

		vList->size = newsize;
		vList->verts = newverts;
	}
	
	for (i=0;i<vList->dim;i++)
	{
		vList->verts[(vList->numVerts)*(vList->dim)+i] = v[i];
	}
	vList->numVerts++;
	return INSERT_SUCCESS;
}

/** @brief Removes a vertex from the list, closes space in memory.
 *
 * @param vList  The vertex list to delete from.
 * @param i Index of the vertex to be deleted.
 * 
 * @retVal DELETE_FAILED Successful insertion.
 * @retVal DELETE_SUCCESS Indicates a failure error code.
 *
 * @note User should be careful using this operation if the vertex list is
 *       index into by an edge list.
 */
r_code DeleteVertex(vertexList *vList,int i)
{
	if (i < 0 || i >= vList->numVerts)
	{
		return DELETE_FAILED;
	}
	else
	{
		int j,s,f;
		void * ptr;
		s = i*(vList->dim);
		f = (vList->numVerts-1)*(vList->dim);
		for (j=s;j<f;j++)
		{
			vList->verts[j] = vList->verts[j + vList->dim];
		}
		ptr = (void *)(vList->verts + f);
		memset(ptr,0,vList->dim*sizeof(float));
		vList->numVerts--;
		return DELETE_SUCCESS;
	}
}

/** @brief Returns the memory address of the ith vertex.
 *
 * @param vList The vertex list.
 * @param i The index of the vertex. 
 * 
 * @returns A pointer to start of the ith vertex.
 *
 * @retVal NULL if it does not exist.
 *
 * @note editing the result of this function will change values in the vertex list
 *       if a copy is desired use GetVertex_cp
 */
float * GetVertex_ptr(vertexList *vList,int i)
{
	if (i < 0 || i >= vList->numVerts)
	{
		return NULL;
	}
	else
	{
		return vList->verts + i*(vList->dim);
	}
}

/** @brief Makes a copy of the ith vertex in the list.
 *
 * @param vList The vertex list.
 * @param i The index of the vertex.
 *
 * @returns A pointer to copy of the ith vertex.
 * @retVal NULL if it does not exist.
 * @retVal OUT_OF_MEMORY If new memory cannot be allocated for the vertex.
 */
float * GetVertex_cp(vertexList *vList,int i)
{
	if (i < 0 || i >= vList->numVerts)
	{
		return NULL;
	}
	else
	{
		float *v;
		int j;
		if(!(v = (float *)malloc(vList->dim*sizeof(float))))
		{
			return OUT_OF_MEMORY;
		}
		
		for(j=0;j<vList->dim;j++)
		{
			v[j] = vList->verts[i*(vList->dim) + j];
		}
		
		return v;
	}
}

/** @brief Allocates memory for face list.
 *
 * @param size The number of faces.
 * @param maxVerts The largest number of vertices allowed per face.
 * 
 * @returns pointer to allocated memory for the face list.
 * @retVal OUT_OF_MEMORY if insufficient memory for creation.
 */
faceList * CreateFaceList(int size,int maxVerts)
{
	faceList *fList;
	int n;
	if (!(fList = (faceList *)malloc(sizeof(faceList))))
	{
		return OUT_OF_MEMORY;
	}
	
	fList->numFaces = 0;
	fList->size = size;
	fList->maxVerts = maxVerts;
	
	n = size*maxVerts;
	
	if (!(fList->faces = (int *)malloc(n*sizeof(int))))
	{
		return OUT_OF_MEMORY;
	}
	
	if (!(fList->faceTypes = (unsigned char *)malloc(size*sizeof(unsigned char))))
	{
		return OUT_OF_MEMORY;
	}
	memset((void *)(fList->faces),0,n*sizeof(int));
	memset((void *)(fList->faceTypes),0,size*sizeof(unsigned char));
	return fList;
}

/** @brief Adds face to the face list.
 *
 * @param fList The face list to insert into.
 * @param face The face to be inserted.
 * @param type The Face element type.
 * 
 * @retVal INSERT_SUCCESS On success.
 * @retVal INSERT_FAILED On failure.
 */
r_code InsertFace(faceList *fList,int *face,unsigned char type)
{
	int j;
	if (fList->size <= fList->numFaces)
	{
		int newsize;
		int *newfaces;
		unsigned char *newtypes;
		newsize = fList->size + REALLOC_EXTRA;
		newfaces = (int*)realloc(fList->faces,newsize*(fList->maxVerts)*sizeof(int));
		if (!newfaces)
		{
			return INSERT_FAILED;
		}
		newtypes = (unsigned char *)realloc(fList->faceTypes,newsize*sizeof(unsigned char));
		if(!newtypes)
		{
			return INSERT_FAILED;
		}

		fList->size = newsize;
		fList->faces = newfaces;
		fList->faceTypes = newtypes;
	}

	for (j=0;j<(int)type;j++)
	{
		fList->faces[(fList->numFaces)*(fList->maxVerts)+j] = face[j];
	}
	fList->faceTypes[fList->numFaces] = type;
	fList->numFaces++;
	return INSERT_SUCCESS;
}

/** @brief Deletes a face from the face list.
 *
 * @param fList The face list to delete from.
 * @param i The index of the face to be deleted.
 * 
 * @retVal DELETE_SUCCESS on success.
 * @retVal DELETE_FAILED  on failure.  
 */
r_code DeleteFace(faceList *fList, int i)
{
	if (i < 0 || i >= fList->numFaces)
	{
		return DELETE_FAILED;
	}
	else
	{
		int j,s,f;
		void * ptr;
		s = i*(fList->maxVerts);
		f = (fList->numFaces -1)*(fList->maxVerts);
		for (j=s;j<f;j++)
		{
			fList->faces[j] = fList->faces[j+fList->maxVerts];
		}
		for (j=i;j<fList->numFaces-1;j++)
		{
			fList->faceTypes[j] = fList->faceTypes[j+1];
		}
		
		ptr = (void *)(fList->faces + f);
		memset(ptr,0,fList->maxVerts*sizeof(float));
		fList->faceTypes[fList->numFaces-1] = 0;
		fList->numFaces--;
		return DELETE_SUCCESS;
	}
}

/** @brief Returns the memory address of the ith face.
 *
 * @param fList The face list.
 * @param i The index of the face. 
 * 
 * @returns The pointer to start of the ith face.
 * @retVal NULL if it does not exist.
 *
 * @note editing the result of this function will change values in the face list
 *       if a copy is desired use GetFace_cp .
 */
int * GetFace_ptr(faceList *fList,int i)
{
	if (i < 0 || i >= fList->numFaces)
	{
		return NULL;
	}
	else
	{
		return fList->faces + i*(fList->maxVerts);
	}
}

/** @brief Returns a copy of the ith face.
 *
 * @param fList The face list.
 * @param i The index of the face. 
 *
 * @returns The pointer to start of the ith face.
 * @retVal NULL if it does not exist.
 * @retVal OUT_OF_MEMORY an error occurred in copy memory allocation.
 *
 */
int * GetFace_cp(faceList *fList,int i)
{
	if (i < 0 || i >= fList->numFaces)
	{
		return NULL;
	}
	else
	{
		int * f;
		int j;
		if (!(f = (int *)malloc(fList->faceTypes[i]*sizeof(int))))
		{
			return OUT_OF_MEMORY;
		}
		for (j=0;j<fList->faceTypes[i];j++)
		{
			f[j] = fList->faces[i*(fList->maxVerts) + j];
		}
		return f;
	}
}

/** @brief Calculates the midpoint of the ith face int the fList.
 *
 * @param fList The face list.
 * @param i The index of the face .
 * @param vList The list of face vertices.
 *
 * @returns A pointer to a the mid-point (geometric centre) of the face.
 * @retVal NULL If an error occurred in calculation.
 * @retVal OUT_OF_MEMORY If an error occurred in memory allocation.
 *
 * 
 */
float * FaceMidPoint(faceList *fList, int i,vertexList *vList)
{
	unsigned char type;
	float *verts;
	float *mean;
	int *face;
	float *vert;
	int j,k;
	
	type = fList->faceTypes[i];
	
	if (!(verts = (float *)malloc(type*(vList->dim)*sizeof(float))))
	{
		return OUT_OF_MEMORY;
	}
	
	face = GetFace_ptr(fList,i);
	
	for (j=0;j<type;j++)
	{
		vert = GetVertex_ptr(vList,face[j]);
		for (k=0;k<(vList->dim);k++)
		{
			verts[j*(vList->dim)+k] = vert[k];
		} 
	}
	
	mean = Mean_f(type,verts,vList->dim);
	
	if (!(mean))
	{
		return NULL;
	}
	
	free(verts);
	
	return mean;
}

/** @brief Calculates the face area of the ith face in the face list.
 *
 * @param fList The face list.
 * @param i The index of face of interest.
 * @param vList The list of vertices indexed by eList.
 *
 * @retruns The area of the face.
 */
float FaceArea(faceList *fList, int i,vertexList *vList)
{
	float *v0,*v1;
	int *face;
	unsigned char type;
	int j;
	float area;
	
	type = fList->faceTypes[i];
	face = GetFace_ptr(fList,i);
	
	area = 0.0;
	for (j=0;j<type-1;j++)
	{
		v0 = GetVertex_ptr(vList,face[j]);
		v1 = GetVertex_ptr(vList,face[j+1]);
		area += MagCross_f(v0,v1);
	}
	
	return area*0.5;
}

/* CreateMesh(): allocates memory for a Polygonal mesh
 *
 * Parameters:
 *     numVerts - number of vertices to allocate for
 *     dim - dimension of vertices
 *     numFaces - number of Faces to allocate for
 *     maxVerts - largest number of vertices any face may have
 * Returns:
 *     m - a pointer to mesh structure, NULL if an error occurred
 */
mesh * CreateMesh(int numVerts,int dim,int numFaces, int maxVerts)
{
	mesh *m;
	
	if (!(m = (mesh *)malloc(sizeof(mesh))))
	{
		return OUT_OF_MEMORY;
	}
	
	m->vList = CreateVertexList(numVerts,dim);
	if (!(m->vList))
	{
		return NULL;
	}
	
	m->fList = CreateFaceList(numFaces,maxVerts);
	if (!(m->fList))
	{
		return NULL;
	}
	return m;
}

/* CopyMesh(): Creates a new mesh structure which is an exact copy for the 
 *             given mesh m
 * Parameter:
 *     m - mesh to copy
 * Returns:
 *    m_cp - mesh identical to m
 */
mesh * CopyMesh(mesh *m)
{
	mesh *m_cp;
	int numVerts,dim,numFaces,maxVerts;
	int i;
	
	numVerts = m->vList->size;
	dim = m->vList->dim;
	numFaces = m->fList->size;
	maxVerts = m->fList->maxVerts;
	
	m_cp = CreateMesh(numVerts,dim,numFaces,maxVerts);
	
	if (!m_cp)
	{
		return NULL;
	}
	
	m_cp->vList->numVerts = m->vList->numVerts;
	m_cp->fList->numFaces = m->fList->numFaces;
	
	memcpy(m_cp->vList->verts,m->vList->verts,numVerts*dim*sizeof(float));
	memcpy(m_cp->fList->faces,m->fList->faces,numFaces*maxVerts*sizeof(int));
	memcpy(m_cp->fList->faceTypes,m->fList->faceTypes,numFaces*sizeof(unsigned char));
	
	return m_cp;
}

/* CreateDual(): Creates the Geometric dual of the given mesh. That is, faces 
 *               become vertices and vertices become faces
 * Parameters:
 *     m - original mesh
 *
 * Returns:
 *    m_dual - geometric dual of m
 */
mesh * CreateDual(mesh *m)
{
	/*TODO: Need some more thought about the most efficient way to do this*/
	mesh *mdual;
	int nv,nf;
	int i,j,k;
	int *face;
	float *face2vert;
	int vert2face[MAX_FACE_TYPE];
	int reorderedvert2face[MAX_FACE_TYPE];
	int edge[2];
	unsigned char type;
	unsigned char newtype;
	unsigned char retype;
	unsigned char in;
	int u,l;
	mdual = CreateMesh(m->fList->numFaces,3,m->vList->numVerts,MAX_FACE_TYPE);
	nf = m->fList->numFaces;
	/*faces become vertices*/
	for (i=0;i<nf;i++)
	{
		face2vert = FaceMidPoint(m->fList,i,m->vList);
		InsertVertex(mdual->vList,face2vert);
		free(face2vert);
	}
	/*vertices become faces*/
	nv = m->vList->numVerts;
	for (j=0;j<nv;j++)
	{
		newtype = 0;
		/*find all faces that use this vertex*/
		for (i=0;i<nf;i++)
		{
			face = GetFace_ptr(m->fList,i);
			type = m->fList->faceTypes[i];
			in = 0;
			/*trying to reduce branchin*/
			for (k=0;k<type;k++)
			{
				in |= (face[k] == j);
			}

			if (in)
			{
				vert2face[newtype] = i;
				newtype++;
			}
		}

		/*TODO: need to apply an vertex reordering here*/
		/*find a face that starts with this vertex*/
		for (i=0;i<newtype;i++)
		{
			face = GetFace_ptr(m->fList,vert2face[i]);
			type = m->fList->faceTypes[vert2face[i]];
			k = 0;
			while (face[k] != j) k++;
					
			u = (k == type-1) ? 0 : k+1;
			edge[0] = j;
			edge[1] = face[u];
			break;
		}
		
		for (retype=0;retype<newtype;retype++)
		{
			for (i=0;i<newtype;i++)
			{
				if (vert2face[i] != -1)
				{
					face = GetFace_ptr(m->fList,vert2face[i]);
					type = m->fList->faceTypes[vert2face[i]];
					k = 0;
					while (face[k] != edge[0]) k++;
					
					u = (k == type-1) ? 0 : k+1;
					l = (k == 0) ? type-1 : k-1;
					if (face[u] == edge[1])
					{
						edge[1] = face[l];
						reorderedvert2face[retype] = vert2face[i];
						vert2face[i] = -1;
						break;
					}
				}
			}
		}

		for (i=0;i<newtype;i++)
		{
			vert2face[i] = reorderedvert2face[i];
		}

		InsertFace(mdual->fList,vert2face,newtype);
	}
	return mdual;
}

/* SubDivideFaces(): subdivides each face by constructing a new face from the
 *                   midpoints of the faces edges
 *
 * Parameters:
 *     m - original mesh
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code SubDivideFaces(mesh *m)
{
	int N,dim;
	int numVerts_prev;
	int i,j;
	int *face;
	float *m_j,*vert_j,*vert_jp1;
	int triface[3];
	r_code rc;
	N = m->fList->numFaces;
	dim = m->vList->dim;
	/*for every face*/
	for (i=0;i<N;i++)
	{
		register int n;
		face = GetFace_ptr(m->fList,i);
		n = m->fList->faceTypes[i];
		numVerts_prev = m->vList->numVerts;
		/*get the midpoints of each edge and insert (indexes are m_j = numVerts_pre + j)*/
		for (j=0;j<n;j++)
		{
			vert_j = GetVertex_ptr(m->vList,face[j]);
			vert_jp1 = GetVertex_ptr(m->vList,face[(j+1)%n]);
			m_j = Midpoint_f(vert_j,vert_jp1,dim);
			rc = InsertVertex(m->vList,m_j);
			if (rc <=0) return rc;
			free(m_j);/*I REALLY HATE THIS!!!! Need to fix the vector math lib*/
		}

		/*insert triangular faces of the form {v_j,m_j,m_(j-1)}*/
		triface[0] = face[0];
		triface[1] = numVerts_prev; 
		triface[2] = numVerts_prev + n-1; 
		rc = InsertFace(m->fList,triface,(unsigned char)3);
		if (rc <=0) return rc;
	
		face = GetFace_ptr(m->fList,i);
		for (j=1;j<n;j++)
		{
			triface[0] = face[j];
			triface[1] = numVerts_prev + j; 
			triface[2] = numVerts_prev + j-1; 
			rc = InsertFace(m->fList,triface,(unsigned char)3);
			if (rc <=0) return rc;
	
			face = GetFace_ptr(m->fList,i);
		}

		/*replace v_j index in face for m_j*/
		for (j=0;j<n;j++)
		{
			face[j] = numVerts_prev + j;
		}
	}
	/*clean the mesh by removing duplicates*/
	rc = RemoveDuplicateVertices(m);
	return rc;
}

/* RemoveDuplicateVertices(): deletes any duplicate vertices int he mesh and 
 *                            fixes any references from faces
 */
r_code RemoveDuplicateVertices(mesh *m)
{
	int *remap;
	int *dups;
	int i,j;
	int N,dim;
	int *numDups;
	int * face;
	float *vert_i,*vert_j;

	N = m->vList->numVerts;
	dim = m->vList->dim;
/*allocate some memory for required arrays*/
	if (!(remap = (int *)malloc(N*sizeof(int))))
	{
		return OUT_OF_MEMORY;
	}
	if (!(dups = (int *)malloc(N*sizeof(int))))
	{
		return OUT_OF_MEMORY;
	}
	if (!(numDups = (int *)malloc(N*sizeof(int))))
	{
		return OUT_OF_MEMORY;
	}
	/*initialise duplicate array to all ones*/
	for (i=0;i<N;i++) dups[i] = -1;
	/*for every vertex*/
	for (i=0;i<N;i++)
	{
		/*if this is not a duplicate*/
		if (dups[i] == -1)
		{
			vert_i = GetVertex_ptr(m->vList,i);
			/*search for duplicates*/
			for (j=i+1;j<N;j++)
			{
				vert_j = GetVertex_ptr(m->vList,j);
				if (NormDiff_f(vert_i,vert_j,dim,TWO_NORM) == 0.0)
				{
					dups[j] = i;
				}
			}
		}
	}
	/*now create remap*/
	numDups[0] = 0;
	for (i=1;i<N;i++)
	{
		numDups[i] = numDups[i-1] + (dups[i-1] != -1);
	}
	for (i=0;i<N;i++)
	{
		/*verts that are not duplicates need to be shifted by the 
		 * number of duplicates so far
		 */
		if (dups[i] == -1)
		{
			remap[i] = i - numDups[i];
		}
		else /*for duplicates assign the original index*/
		{
			remap[i] = dups[i]-numDups[dups[i]];
			/*delete the duplicate*/
			DeleteVertex(m->vList,i - numDups[i]);
		}
	}
	/*re-map faces to new index values*/
	N = (m->fList->numFaces);
	for (i=0;i<N;i++)
	{
		face = GetFace_ptr(m->fList,i);
		for (j=0;j<m->fList->faceTypes[i];j++)
		{
			face[j] = remap[face[j]];
		}
	}
	/*clean up*/
	free(remap);
	free(dups);	
	free(numDups);
	return SUCCESS; 
}


/* SurfaceArea(): Calculates the surface area of the given mesh. Simply the sum 
 *                the Face surface areas
 * Parameters:
 *      m - the mesh to compute the are of
 */
float SurfaceArea(mesh* m)
{
	float surfArea;
	int i;
	
	surfArea = 0.0;
	for (i=0;i<m->fList->numFaces;i++)
	{
		surfArea += FaceArea(m->fList,i,m->vList);
	}
	
	return surfArea;
}

/* GeometricCentre(): Calculates the geometric centre of the mesh
 *
 * Parameters:
 *      m - mesh to get centre of
 */
float * GeometricCentre(mesh* m)
{
	float * geoCentre;
	
	geoCentre = Mean_f(m->vList->numVerts,m->vList->verts,m->vList->dim);
	
	return geoCentre;
}

/* CreateIcosahedron(): does what is says it does...
 *
 * Parameters:
 *      none...
 * Returns:
 *    a mesh that resembles an icosahedron
 */
mesh * CreateIcosahedron(void)
{
	mesh *m;
	float vert[3];
	int face[3];
	/*12 verts and 20 faces*/
	m = CreateMesh(12,3,20,3);
	
	vert[0]=0.0;vert[1]=0.0;vert[2]=2.0;
	InsertVertex(m->vList,vert);
	vert[0]=1.788854;vert[1]= 0.000000;vert[2]= 0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=0.552786;	vert[1]= 1.701302;	vert[2]= 0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=-1.447214;	vert[1]= 1.051462;	vert[2]= 0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=-1.447214;	vert[1]= -1.051462;	vert[2]= 0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=0.552786;	vert[1]= -1.701302;	vert[2]= 0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=1.447214;	vert[1]= 1.051462;	vert[2]= -0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=-0.552786;	vert[1]= 1.701302;	vert[2]= -0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=-1.788854;	vert[1]= 0.000000;	vert[2]= -0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=-0.552786;	vert[1]= -1.701302;	vert[2]= -0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=1.447214;	vert[1]= -1.051462;	vert[2]= -0.894427;
	InsertVertex(m->vList,vert);
	vert[0]=0.0;	vert[1]= 0.0;	vert[2]= -2.0;
	InsertVertex(m->vList,vert);
	face[0]= 2;	face[1]= 0;	face[2]= 1;
	InsertFace(m->fList,face,3);
	face[0]= 3;	face[1]= 0;	face[2]= 2;
	InsertFace(m->fList,face,3);
	face[0]= 4;	face[1]= 0;	face[2]= 3;
	InsertFace(m->fList,face,3);
	face[0]= 5;	face[1]= 0;	face[2]= 4;
	InsertFace(m->fList,face,3);
	face[0]= 1;	face[1]= 0;	face[2]= 5;
	InsertFace(m->fList,face,3);
	face[0]= 2;	face[1]= 1;	face[2]= 6;
	InsertFace(m->fList,face,3);
	face[0]= 7;	face[1]= 2;	face[2]= 6;
	InsertFace(m->fList,face,3);
	face[0]= 3;	face[1]= 2;	face[2]= 7;
	InsertFace(m->fList,face,3);
	face[0]= 8;	face[1]= 3;	face[2]= 7;
	InsertFace(m->fList,face,3);
	face[0]= 4;	face[1]= 3;	face[2]= 8;
	InsertFace(m->fList,face,3);
	face[0]= 9;	face[1]= 4;	face[2]= 8;
	InsertFace(m->fList,face,3);
	face[0]= 5;	face[1]= 4;	face[2]= 9;
	InsertFace(m->fList,face,3);
	face[0]= 10;	face[1]= 5;	face[2]= 9;
	InsertFace(m->fList,face,3);
	face[0]= 6;	face[1]= 1 ;	face[2]=10;
	InsertFace(m->fList,face,3);
	face[0]= 1;	face[1]= 5;	face[2]= 10;
	InsertFace(m->fList,face,3);
	face[0]= 6;	face[1]= 11;	face[2]= 7;
	InsertFace(m->fList,face,3);
	face[0]= 7;	face[1]= 11;	face[2]= 8;
	InsertFace(m->fList,face,3);
	face[0]= 8;	face[1]= 11;	face[2]= 9;
	InsertFace(m->fList,face,3);
	face[0]= 9;	face[1]= 11;	face[2]= 10;
	InsertFace(m->fList,face,3);
	face[0]= 10 ;	face[1]=11;	face[2]= 6;
	InsertFace(m->fList,face,3);

	return m;
}

/* CreateTorus(): does what is says it does...
 *
 * Parameters:
 *      none...
 * Returns:
 *    a mesh that resembles an torus
 */
mesh * CreateTorus(void)
{
	mesh *m;
	float vert[3];
	int face[3];
	/*40 verts and 20 faces*/
	m = CreateMesh(40,3,80,3);

	vert[0] = 1.000000;	vert[1] = 0.000000;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.309000;	vert[1] = 0.951100;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.809000;	vert[1] = 0.587800;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.809000;	vert[1] = -0.587800;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.309000;	vert[1] = -0.951100;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.401300;	vert[1] = 1.018100;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.535200;	vert[1] = 1.647300;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.732100;	vert[1] = 0.000000;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.535200;	vert[1] = -1.647300;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.401300;	vert[1] = -1.018100;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.200650;	vert[1] = 0.509050;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.855150;	vert[1] = 0.984600;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.654500;	vert[1] = 0.475550;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.113100;	vert[1] = 1.299200;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.672100;	vert[1] = 1.117550;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.250000;	vert[1] = 0.769450;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.270550;	vert[1] = 0.293900;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.270550;	vert[1] = -0.293900;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.809000;	vert[1] = 0.000000;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.672100;	vert[1] = -1.117550;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.113100;	vert[1] = -1.299200;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.250000;	vert[1] = -0.769450;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.855150;	vert[1] = -0.984600;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.200650;	vert[1] = -0.509050;	vert[2] = 0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.654500;	vert[1] = -0.475550;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.433050;	vert[1] = 1.332700;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.133650;	vert[1] = 0.823650;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.133650;	vert[1] = -0.823650;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.433050;	vert[1] = -1.332700;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.401300;	vert[1] = 0.000000;	vert[2] = 0.000000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.200650;	vert[1] = 0.509050;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.855150;	vert[1] = 0.984600;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.113100;	vert[1] = 1.299200;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.672100;	vert[1] = 1.117550;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.270550;	vert[1] = 0.293900;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -1.270550;	vert[1] = -0.293900;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.672100;	vert[1] = -1.117550;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = -0.113100;	vert[1] = -1.299200;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 0.855150;	vert[1] = -0.984600;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);
	vert[0] = 1.200650;	vert[1] = -0.509050;	vert[2] = -0.200000;
	InsertVertex(m->vList,vert);

	face[0]=10;	face[1]=11;	face[2]=12;
	InsertFace(m->fList,face,3);
	face[0]=13;	face[1]=14;	face[2]=15;
	InsertFace(m->fList,face,3);
	face[0]=16;	face[1]=17;	face[2]=18;
	InsertFace(m->fList,face,3);
	face[0]=19;	face[1]=20;	face[2]=21;
	InsertFace(m->fList,face,3);
	face[0]=22;	face[1]=23;	face[2]=24;
	InsertFace(m->fList,face,3);
	face[0]=25;	face[1]=13;	face[2]=11;
	InsertFace(m->fList,face,3);
	face[0]=26;	face[1]=16;	face[2]=14;
	InsertFace(m->fList,face,3);
	face[0]=27;	face[1]=19;	face[2]=17;
	InsertFace(m->fList,face,3);
	face[0]=28;	face[1]=22;	face[2]=20;
	InsertFace(m->fList,face,3);
	face[0]=29;	face[1]=10;	face[2]=23;
	InsertFace(m->fList,face,3);
	face[0]=0;	face[1]=10;	face[2]=12;
	InsertFace(m->fList,face,3);
	face[0]=5;	face[1]=11;	face[2]=10;
	InsertFace(m->fList,face,3);
	face[0]=1;	face[1]=12;	face[2]=11;
	InsertFace(m->fList,face,3);
	face[0]=1;	face[1]=13;	face[2]=15;
	InsertFace(m->fList,face,3);
	face[0]=6;	face[1]=14;	face[2]=13;
	InsertFace(m->fList,face,3);
	face[0]=2;	face[1]=15;	face[2]=14;
	InsertFace(m->fList,face,3);
	face[0]=2;	face[1]=16;	face[2]=18;
	InsertFace(m->fList,face,3);
	face[0]=7;	face[1]=17;	face[2]=16;
	InsertFace(m->fList,face,3);
	face[0]=3;	face[1]=18;	face[2]=17;
	InsertFace(m->fList,face,3);
	face[0]=3;	face[1]=19;	face[2]=21;
	InsertFace(m->fList,face,3);
	face[0]=8;	face[1]=20;	face[2]=19;
	InsertFace(m->fList,face,3);
	face[0]=4;	face[1]=21;	face[2]=20;
	InsertFace(m->fList,face,3);
	face[0]=4;	face[1]=22;	face[2]=24;
	InsertFace(m->fList,face,3);
	face[0]=9;	face[1]=23;	face[2]=22;
	InsertFace(m->fList,face,3);
	face[0]=0;	face[1]=24;	face[2]=23;
	InsertFace(m->fList,face,3);
	face[0]=5;	face[1]=25;	face[2]=11;
	InsertFace(m->fList,face,3);
	face[0]=6;	face[1]=13;	face[2]=25;
	InsertFace(m->fList,face,3);
	face[0]=1;	face[1]=11;	face[2]=13;
	InsertFace(m->fList,face,3);
	face[0]=6;	face[1]=26;	face[2]=14;
	InsertFace(m->fList,face,3);
	face[0]=7;	face[1]=16;	face[2]=26;
	InsertFace(m->fList,face,3);
	face[0]=2;	face[1]=14;	face[2]=16;
	InsertFace(m->fList,face,3);
	face[0]=7;	face[1]=27;	face[2]=17;
	InsertFace(m->fList,face,3);
	face[0]=8;	face[1]=19;	face[2]=27;
	InsertFace(m->fList,face,3);
	face[0]=3;	face[1]=17;	face[2]=19;
	InsertFace(m->fList,face,3);
	face[0]=8;	face[1]=28;	face[2]=20;
	InsertFace(m->fList,face,3);
	face[0]=9;	face[1]=22;	face[2]=28;
	InsertFace(m->fList,face,3);
	face[0]=4;	face[1]=20;	face[2]=22;
	InsertFace(m->fList,face,3);
	face[0]=9;	face[1]=29;	face[2]=23;
	InsertFace(m->fList,face,3);
	face[0]=5;	face[1]=10;	face[2]=29;
	InsertFace(m->fList,face,3);
	face[0]=0;	face[1]=23;	face[2]=10;
	InsertFace(m->fList,face,3);
	face[0]=30;	face[1]=12;	face[2]=31;
	InsertFace(m->fList,face,3);
	face[0]=32;	face[1]=15;	face[2]=33;
	InsertFace(m->fList,face,3);
	face[0]=34;	face[1]=18;	face[2]=35;
	InsertFace(m->fList,face,3);
	face[0]=36;	face[1]=21;	face[2]=37;
	InsertFace(m->fList,face,3);
	face[0]=38;	face[1]=24;	face[2]=39;
	InsertFace(m->fList,face,3);
	face[0]=25;	face[1]=31;	face[2]=32;
	InsertFace(m->fList,face,3);
	face[0]=26;	face[1]=33;	face[2]=34;
	InsertFace(m->fList,face,3);
	face[0]=27;	face[1]=35;	face[2]=36;
	InsertFace(m->fList,face,3);
	face[0]=28;	face[1]=37;	face[2]=38;
	InsertFace(m->fList,face,3);
	face[0]=29;	face[1]=39;	face[2]=30;
	InsertFace(m->fList,face,3);
	face[0]=30;	face[1]=0;	face[2]=12;
	InsertFace(m->fList,face,3);
	face[0]=31;	face[1]=5;	face[2]=30;
	InsertFace(m->fList,face,3);
	face[0]=12;	face[1]=1;	face[2]=31;
	InsertFace(m->fList,face,3);
	face[0]=32;	face[1]=1;	face[2]=15;
	InsertFace(m->fList,face,3);
	face[0]=33;	face[1]=6;	face[2]=32;
	InsertFace(m->fList,face,3);
	face[0]=15;	face[1]=2;	face[2]=33;
	InsertFace(m->fList,face,3);
	face[0]=34;	face[1]=2;	face[2]=18;
	InsertFace(m->fList,face,3);
	face[0]=35;	face[1]=7;	face[2]=34;
	InsertFace(m->fList,face,3);
	face[0]=18;	face[1]=3;	face[2]=35;
	InsertFace(m->fList,face,3);
	face[0]=36;	face[1]=3;	face[2]=21;
	InsertFace(m->fList,face,3);
	face[0]=37;	face[1]=8;	face[2]=36;
	InsertFace(m->fList,face,3);
	face[0]=21;	face[1]=4;	face[2]=37;
	InsertFace(m->fList,face,3);
	face[0]=38;	face[1]=4;	face[2]=24;
	InsertFace(m->fList,face,3);
	face[0]=39;	face[1]=9;	face[2]=38;
	InsertFace(m->fList,face,3);
	face[0]=24;	face[1]=0;	face[2]=39;
	InsertFace(m->fList,face,3);
	face[0]=25;	face[1]=5;	face[2]=31;
	InsertFace(m->fList,face,3);
	face[0]=32;	face[1]=6;	face[2]=25;
	InsertFace(m->fList,face,3);
	face[0]=31;	face[1]=1;	face[2]=32;
	InsertFace(m->fList,face,3);
	face[0]=26;	face[1]=6;	face[2]=33;
	InsertFace(m->fList,face,3);
	face[0]=34;	face[1]=7;	face[2]=26;
	InsertFace(m->fList,face,3);
	face[0]=33;	face[1]=2;	face[2]=34;
	InsertFace(m->fList,face,3);
	face[0]=27;	face[1]=7;	face[2]=35;
	InsertFace(m->fList,face,3);
	face[0]=36;	face[1]=8;	face[2]=27;
	InsertFace(m->fList,face,3);
	face[0]=35;	face[1]=3;	face[2]=36;
	InsertFace(m->fList,face,3);
	face[0]=28;	face[1]=8;	face[2]=37;
	InsertFace(m->fList,face,3);
	face[0]=38;	face[1]=9;	face[2]=28;
	InsertFace(m->fList,face,3);
	face[0]=37;	face[1]=4;	face[2]=38;
	InsertFace(m->fList,face,3);
	face[0]=29;	face[1]=9;	face[2]=39;
	InsertFace(m->fList,face,3);
	face[0]=30; face[1]=5;	face[2]=29;
	InsertFace(m->fList,face,3);
	face[0]=39;	face[1]=0;	face[2]=30;
	InsertFace(m->fList,face,3);
	return m;	
}

/* CreateDoubleTorus(): does what is says it does...
 *
 * Parameters:
 *      none...
 * Returns:
 *    a mesh that resembles a double torus
 */
mesh * CreateDoubleTorus(void)
{
	mesh *m;
	float vert[3];
	int face[3];
	m = CreateMesh(40,3,84,3);
	vert[0]=-0.418007;vert[1]= -0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.306002;vert[1]= 0.000001;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.418006;vert[1]= 0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.642016;vert[1]= 0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.754021;vert[1]= -0.000001;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.642017;vert[1]= -0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.000000;vert[1]= -0.306003;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.000000;vert[1]= 0.306002;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.530012;vert[1]= 0.612005;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.060023;vert[1]= 0.306003;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.060023;vert[1]= -0.306002;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.530011;vert[1]= -0.612005;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.030012;vert[1]= 0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.030012;vert[1]= 0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.530011;vert[1]= 0.000000;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.030011;vert[1]= -0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.030012;vert[1]= -0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.418007;vert[1]= -0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.306002;vert[1]= 0.000001;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.418006;vert[1]= 0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.642016;vert[1]= 0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.754021;vert[1]= -0.000001;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.642017;vert[1]= -0.193997;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.530012;vert[1]= 0.612005;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=1.060023;vert[1]= 0.306003;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=1.060023;vert[1]= -0.306002;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=0.530011;vert[1]= -0.612005;vert[2]= 0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=1.030012;vert[1]= 0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=1.530011;vert[1]= 0.000000;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=1.030011;vert[1]= -0.866025;vert[2]= 0.000000;
	InsertVertex(m->vList,vert);
	vert[0]=0.000000;vert[1]= -0.306003;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.000000;vert[1]= 0.306002;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.530012;vert[1]= 0.612005;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.060023;vert[1]= 0.306003;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-1.060023;vert[1]= -0.306002;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=-0.530011;vert[1]= -0.612005;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=0.530012;vert[1]= 0.612005;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=1.060023;vert[1]= 0.306003;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=1.060023;vert[1]= -0.306002;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	vert[0]=0.530011;vert[1]= -0.612005;vert[2]= -0.200000;
	InsertVertex(m->vList,vert);
	 face[0]=0;face[1]= 6;face[2]= 1;
	InsertFace(m->fList,face,3);
	 face[0]=1;face[1]= 6;face[2]= 7;
	InsertFace(m->fList,face,3);
	 face[0]=1;face[1]= 7;face[2]= 2;
	InsertFace(m->fList,face,3);
	 face[0]=2;face[1]= 7;face[2]= 8;
	InsertFace(m->fList,face,3);
	 face[0]=2;face[1]= 8;face[2]= 3;
	InsertFace(m->fList,face,3);
 	face[0]=3;face[1]= 8;face[2]= 9;
	InsertFace(m->fList,face,3);
	 face[0]=3;face[1]= 9;face[2]= 4;
	InsertFace(m->fList,face,3);
	 face[0]=4;face[1]= 9;face[2]= 10 ;
	InsertFace(m->fList,face,3);
	 face[0]=4;face[1]= 10;face[2]= 5;
	InsertFace(m->fList,face,3);
	 face[0]=5;face[1]= 10;face[2]= 11;
	InsertFace(m->fList,face,3);
	 face[0]=5;face[1]= 11;face[2]= 0;
	InsertFace(m->fList,face,3);
	 face[0]=0;face[1]= 11;face[2]= 6;
	InsertFace(m->fList,face,3);
	 face[0]=7;face[1]= 12;face[2]= 8;
	InsertFace(m->fList,face,3);
	 face[0]=8;face[1]= 12;face[2]= 13;
	InsertFace(m->fList,face,3);
	 face[0]=8;face[1]= 13;face[2]= 9;
	InsertFace(m->fList,face,3);
	 face[0]=9;face[1]= 13;face[2]= 14;
	InsertFace(m->fList,face,3);
	 face[0]=9;face[1]= 14;face[2]= 10;
	InsertFace(m->fList,face,3);
	 face[0]=10;face[1]= 14;face[2]= 15;
	InsertFace(m->fList,face,3);
	 face[0]=10;face[1]= 15;face[2]= 11;
	InsertFace(m->fList,face,3);
	 face[0]=11;face[1]= 15;face[2]= 16;
	InsertFace(m->fList,face,3);
	 face[0]= 11;face[1]= 16;face[2]= 6;
	InsertFace(m->fList,face,3);
	 face[0]=17;face[1]= 6;face[2]= 26;
	InsertFace(m->fList,face,3);
	 face[0]=22;face[1]= 17;face[2]= 26;
	InsertFace(m->fList,face,3);
	 face[0]=22;face[1]= 26;face[2]= 25;
	InsertFace(m->fList,face,3);
	 face[0]=21;face[1]= 22;face[2]= 25;
	InsertFace(m->fList,face,3);
	 face[0]=21;face[1]= 25;face[2]= 24;
	InsertFace(m->fList,face,3);
	 face[0]=20;face[1]= 21;face[2]= 24;
	InsertFace(m->fList,face,3);
	 face[0]=20;face[1]= 24;face[2]= 23;
	InsertFace(m->fList,face,3);
	 face[0]=19;face[1]= 20;face[2]= 23;
	InsertFace(m->fList,face,3);
	 face[0]=19;face[1]= 23;face[2]= 7;
	InsertFace(m->fList,face,3);
	 face[0]=18;face[1]= 19;face[2]= 7;
	InsertFace(m->fList,face,3);
	 face[0]=18;face[1]= 7;face[2]= 6;
	InsertFace(m->fList,face,3);
	 face[0]=17;face[1]= 18;face[2]= 6;
	InsertFace(m->fList,face,3);
	 face[0]=26;face[1]= 6;face[2]= 16;
	InsertFace(m->fList,face,3);
	 face[0]=26;face[1]= 16;face[2]= 29;
	InsertFace(m->fList,face,3);
	 face[0]=25;face[1]= 26;face[2]= 29;
	InsertFace(m->fList,face,3);
	 face[0]=25;face[1]= 29;face[2]= 28;
	InsertFace(m->fList,face,3);
	 face[0]=24;face[1]= 25;face[2]= 28;
	InsertFace(m->fList,face,3);
	 face[0]=24;face[1]= 28;face[2]= 27;
	InsertFace(m->fList,face,3);
	 face[0]=23;face[1]= 24;face[2]= 27;
	InsertFace(m->fList,face,3);
	 face[0]=23;face[1]= 27;face[2]= 12;
	InsertFace(m->fList,face,3);
	 face[0]=7;face[1]= 23;face[2]= 12;
	InsertFace(m->fList,face,3);
	 face[0]=0;face[1]= 1;face[2]= 30;
	InsertFace(m->fList,face,3);
	 face[0]=1;face[1]= 31;face[2]= 30;
	InsertFace(m->fList,face,3);
	 face[0]=1;face[1]= 2;face[2]= 31;
	InsertFace(m->fList,face,3);
	 face[0]=2;face[1]= 32;face[2]= 31;
	InsertFace(m->fList,face,3);
	 face[0]=2;face[1]= 3;face[2]= 32;
	InsertFace(m->fList,face,3);
	 face[0]=3;face[1]= 33;face[2]= 32;
	InsertFace(m->fList,face,3);
	 face[0]=3;face[1]= 4;face[2]= 33;
	InsertFace(m->fList,face,3);
	 face[0]=4;face[1]= 34;face[2]= 33;
	InsertFace(m->fList,face,3);
	 face[0]=4;face[1]= 5 ;face[2]=34;
	InsertFace(m->fList,face,3);
	 face[0]= 5;face[1]= 35;face[2]= 34;
	InsertFace(m->fList,face,3);
	 face[0]=5;face[1]= 0;face[2]= 35;
	InsertFace(m->fList,face,3);
	 face[0]=0;face[1]= 30;face[2]= 35;
	InsertFace(m->fList,face,3);
	 face[0]=31;face[1]= 32;face[2]= 12;
	InsertFace(m->fList,face,3);
	 face[0]=32;face[1]= 13;face[2]= 12;
	InsertFace(m->fList,face,3);
	 face[0]=32;face[1]= 33;face[2]= 13;
	InsertFace(m->fList,face,3);
	 face[0]=33;face[1]= 14;face[2]= 13;
	InsertFace(m->fList,face,3);
	 face[0]=33;face[1]= 34;face[2]= 14;
	InsertFace(m->fList,face,3);
	 face[0]=34;face[1]= 15;face[2]= 14;
	InsertFace(m->fList,face,3);
	 face[0]=34;face[1]= 35;face[2]= 15;
	InsertFace(m->fList,face,3);
	 face[0]=35 ;face[1]=16;face[2]= 15;
	InsertFace(m->fList,face,3);
	 face[0]=35;face[1]= 30;face[2]= 16;
	InsertFace(m->fList,face,3);
	 face[0]=17;face[1]= 39;face[2]= 30;
	InsertFace(m->fList,face,3);
	 face[0]=22;face[1]= 39;face[2]= 17;
	InsertFace(m->fList,face,3);
	 face[0]=22;face[1]= 38;face[2]= 39;
	InsertFace(m->fList,face,3);
	 face[0]=21;face[1]= 38;face[2]= 22;
	InsertFace(m->fList,face,3);
	 face[0]=21;face[1]= 37;face[2]= 38;
	InsertFace(m->fList,face,3);
	 face[0]= 20;face[1]= 37;face[2]= 21;
	InsertFace(m->fList,face,3);
	 face[0]=20;face[1]= 36;face[2]= 37;
	InsertFace(m->fList,face,3);
	 face[0]=19;face[1]= 36;face[2]= 20;
	InsertFace(m->fList,face,3);
	 face[0]=19;face[1]= 31;face[2]= 36;
	InsertFace(m->fList,face,3);
	 face[0]=18;face[1]= 31;face[2]= 19;
	InsertFace(m->fList,face,3);
	 face[0]=18;face[1]= 30;face[2]= 31;
	InsertFace(m->fList,face,3);
	 face[0]=17;face[1]= 30;face[2]= 18;
	InsertFace(m->fList,face,3);
	 face[0]=39;face[1]= 16;face[2]= 30;
	InsertFace(m->fList,face,3);
	 face[0]=39;face[1]= 29;face[2]= 16;
	InsertFace(m->fList,face,3);
	 face[0]=38;face[1]= 29;face[2]= 39;
	InsertFace(m->fList,face,3);
	 face[0]=38;face[1]= 28;face[2]= 29;
	InsertFace(m->fList,face,3);
	 face[0]=37;face[1]= 28;face[2]= 38;
	InsertFace(m->fList,face,3);
	 face[0]=37;face[1]= 27;face[2]= 28;
	InsertFace(m->fList,face,3);
	 face[0]=36;face[1]= 27;face[2]= 37;
	InsertFace(m->fList,face,3);
	 face[0]=36;face[1]= 12;face[2]= 27;
	InsertFace(m->fList,face,3);
	 face[0]=31;face[1]= 12;face[2]= 36;
	InsertFace(m->fList,face,3);
	return m;
}

/* CreateMeshTopology(): Creates a mesh that with at most numFaces Faces, 
 *                   that is of the algebraic topology class if the 
 *                   given genus.
 * Parameters:
 *      numFaces - the minimum number of faces
 *      genus - the genus of the topology (i.e., number of holes)
 * Returns:
 *      a mesh of the given genus
 * NOTE: For now, only genus-0, genus-1, genus-2 are supported
 * TODO: generalise to genus-n
 */
mesh * CreateMeshTopology(int numFaces,int genus)
{
	mesh *m;
	switch(genus)
	{
		case 0: /*sphere*/
			/*start with an icosahedron*/
			m = CreateIcosahedron();
			while(m->fList->numFaces < numFaces)
			{
				SubDivideFaces(m);			
			}
			break;
		case 1:
			m = CreateTorus();
			while(m->fList->numFaces < numFaces)
			{
				SubDivideFaces(m);			
			}
			break;
		case 2:
			m = CreateDoubleTorus();
			while(m->fList->numFaces < numFaces)
			{
				SubDivideFaces(m);			
			}
			break;
		default:
			m = NULL;
			break;
	}
	return m;
}

/* LoadMesh(): reads a mesh from given format
 *
 * Parameters: 
 *     filename - the name of input file
 *     format - the format of the input file
 *
 * Returns:
 *     pointer to a valid mesh if successful, otherwise NULL.
 */
mesh * LoadMesh(char *filename,unsigned char format)
{
	mesh *m;
	int rc;
	switch (format)
	{
		default:
		case OFF_FORMAT:
			rc = ReadOFF(filename,&m);
			break;
		case OBJ_FORMAT:
			rc = ReadOBJ(filename,&m);
			break;
		case STL_FORMAT:
			rc = ReadSTL(filename,&m);
			break;
		case VRML_FORMAT:
			rc = ReadVRML(filename,&m);
			break;
	}
	if (CheckErr(rc))
	{
		return NULL;
	}
	return m;
}

/* SaveMesh(): writes mesh to given format
 *
 * Parameters: 
 *     filename - the name of output file
 *     format - the format of the output file
 *
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code SaveMesh(char *filename, mesh *m,unsigned char format)
{
	switch (format)
	{
		default:
		case OFF_FORMAT:
			return WriteOFF(filename,m);
			break;
		case OBJ_FORMAT:
			return WriteOBJ(filename,m);
			break;
		case STL_FORMAT:
			return WriteSTL(filename,m);
			break;
		case VRML_FORMAT:
			return WriteVRML(filename,m);
			break;
	}
}

/* WriteOFF(): writes mesh in *.off format
 *
 * Parameters:
 *     filename -  then name of the .off file
 *     m - the mesh structure to write to file
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code WriteOFF(char *filename,mesh *m)
{

	FILE* fp;
	int i,j,nv,nf,ne;
	float *vert;
	int *face;
	int lastvert;
	nv = m->vList->numVerts;
	nf = m->fList->numFaces;

	ne = 0;
	for (i=0;i<nf;i++)
	{
		ne += (int)(m->fList->faceTypes[i]);
	}

	if (!(fp = fopen(filename,"w")))
	{
		return WRITE_FAILED;
	}
	fprintf(fp,"OFF\n");
	fprintf(fp,"# no comment :)\n");
	fprintf(fp,"%d %d %d\n",nv,nf,ne);
	for (i=0;i<nv;i++)
	{
		vert = GetVertex_ptr(m->vList,i); 
		fprintf(fp,"%f %f %f\n",vert[X],vert[Y],vert[Z]);
	}
	for (i=0;i<nf;i++)
	{
		face = GetFace_ptr(m->fList,i);
		fprintf(fp,"%u",m->fList->faceTypes[i]);

		for (j=0;j<m->fList->faceTypes[i];j++)
		{
			fprintf(fp," %d",face[j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	return WRITE_SUCCESS;
}

/* WriteOBJ(): writes mesh in *.obj format
 *
 * Parameters:
 *     filename -  then name of the .obj file
 *     m - the mesh structure to write to file
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code WriteOBJ(char *filename,mesh *m)
{
	FILE* fp;
	int i,j,nv,nf,ne;
	float *vert;
	int *face;

	nv = m->vList->numVerts;
	nf = m->fList->numFaces;
	ne = 0;
	for (i=0;i<nf;i++)
	{
		ne += (int)(m->fList->faceTypes[i]);
	}

	if (!(fp = fopen(filename,"w")))
	{
		return WRITE_FAILED;
	}
	fprintf(fp,"# no comment :)\n");
	for (i=0;i<nv;i++)
	{
		vert = GetVertex_ptr(m->vList,i); 
		fprintf(fp,"v %f %f %f\n",vert[X],vert[Y],vert[Z]);
	}
	for (i=0;i<nf;i++)
	{
		face = GetFace_ptr(m->fList,i);
		fprintf(fp,"f");
		for (j=0;j<m->fList->faceTypes[i];j++)
		{
			fprintf(fp," %d",face[j]+1);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	return WRITE_SUCCESS;
}

/* WriteSTL(): writes mesh in *.stl format
 *
 * Parameters:
 *     filename -  then name of the .stl file
 *     m - the mesh structure to write to file
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code WriteSTL(char *filename,mesh *m)
{
	FILE* fp;
	int i,j,nv,nf,ne;
	float *vert;
	int *face;
	float *v0,*v1,*vn;
	int * e0,*e1,*en;
	float *norm;
	
	if (!(fp = fopen(filename,"w")))
	{
		return WRITE_FAILED;
	}
	
	nv = m->vList->numVerts;
	nf = m->fList->numFaces;
	ne = 0;
	for (i=0;i<nf;i++)
	{
		ne += (int)(m->fList->faceTypes[i]);
	}
	
	fprintf(fp,"solid mesh\n");
	for (i=0;i<nf;i++)
	{
		face = GetFace_ptr(m->fList,i);
		v0 = GetVertex_ptr(m->vList,face[0]);
		v1 = GetVertex_ptr(m->vList,face[1]);
		vn = GetVertex_ptr(m->vList,face[m->fList->faceTypes[i]-1]);
				
		norm = Normal_f(v0, v1,vn);
		fprintf(fp,"facet ");
		fprintf(fp,"normal %f %f %f\n",norm[X],norm[Y],norm[Z]);
		fprintf(fp,"outer loop\n");
		for (j=0;j<m->fList->faceTypes[i];j++)
		{
			vert = GetVertex_ptr(m->vList,face[j]);
			fprintf(fp,"vertex %f %f %f\n",vert[X],vert[Y],vert[Z]);
		}
		fprintf(fp,"endloop\n");
		fprintf(fp,"endfacet\n");
	}
	fprintf(fp,"endsolid mesh");
	fclose(fp);

	return WRITE_SUCCESS;
}

/* WriteVRML(): writes mesh in *.vrml format
 *
 * Parameters:
 *     filename -  then name of the .vrml file
 *     m - the mesh structure to write to file
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate error code.
 */
r_code WriteVRML(char *filename,mesh *m)
{
	return NOT_IMPLEMENTED;
}

/* ReadOFF(): imports a mesh from a .off file... 
 *
 * Parameters:
 *     filename - I hope you are intelligent enought to work
 *                this one out
 *     m - a pointer to store the address of the new mesh
 * Returns:
 *     r_code - > 0 on success and <= 0 on fail with appropriate
 *     error code. Like most functions in this library...
 */
r_code ReadOFF(char *filename, mesh **m)
{
	FILE *fp;
	int nv,nf,ne;
	float vert[3];
	int face[8];
	char c;
	char buffer[25];
	int i,j,k;
	int faceVerts;
	int maxVerts;
    fpos_t pos;	
	if (!(fp = fopen(filename,"r")))
	{
		return READ_FAILED; 	
	}
	/*first like must be this*/
	fscanf(fp,"%s\n",buffer);
	if ((buffer[0] != 'O') || (buffer[1] != 'F') || (buffer[2] != 'F'))
	{
		return READ_FAILED;
	}

	c = fgetc(fp);
	while (c == '#')
	{
		/*read to the end of the line*/
		while(c != '\n') c = fgetc(fp);
		/*check next line*/
		c = fgetc(fp);
	}
	/*rewind one step*/
	fseek(fp,-1,SEEK_CUR);
	/*read numVerts numFaces and numEdges*/
	fscanf(fp,"%d %d %d\n",&nv,&nf,&ne);
	/*find the max number of vertices per face*/
	fgetpos(fp,&pos);
	for(i=0;i<nv;i++)
	{
		c = fgetc(fp);
		while(c != '\n') c = fgetc(fp);
	}
	maxVerts = 0;
	for (i=0;i<nf;i++)
	{
		fscanf(fp,"%d",&faceVerts);
		if (maxVerts < faceVerts) maxVerts = faceVerts;
		c = fgetc(fp);
		while(c != '\n') c = fgetc(fp);

	}
	fsetpos(fp,&pos);/*return the beginning of data*/
	*m = CreateMesh(nv,3,nf,maxVerts);
	/*read and store vertices*/
	for (i=0;i<nv;i++)
	{
		fscanf(fp,"%f %f %f\n",vert,vert+1,vert+2);
		InsertVertex((*m)->vList,vert);
	}

	/*now read faces*/
	for (i=0;i<nf;i++)
	{
		/*n,v0,v1,v2,...,v(n-1)[,R,G,B,A]*/
		fscanf(fp,"%d",&faceVerts);
		for (j=0;j<faceVerts;j++)
		{
			fscanf(fp,"%d",face+j);
		}
		/*re to end of the line... in case of RGBA data*/
		c = fgetc(fp);
		while(c != '\n') c = fgetc(fp);
		/*insert the face*/
		InsertFace((*m)->fList,face,(unsigned char)faceVerts);
	}
	fclose(fp);	
	return READ_SUCCESS;
}


r_code ReadOBJ(char *filename,mesh **m)
{
	return NOT_IMPLEMENTED;
}

r_code ReadSTL(char *filename,mesh **m)
{
	return NOT_IMPLEMENTED;
}

r_code ReadVRML(char *filename, mesh **m)
{
	return NOT_IMPLEMENTED;
}

/* CheckErr(): Checks if the return code contains an error code
 *
 * Parameters: 
 *     rc - return code to test
 * Returns:
 *     1 if error has occurred, 0 otherwise
 */
unsigned char CheckErr(r_code rc)
{
	if (rc <= 0 )
	{
		PrintErrorMsg(rc);
		return 1;
	}
	return 0;
}

/* PrintErrorMsg(): Prints error message to stderr
 *
 * Parameters:
 *     err - error code
 */
void PrintErrorMsg(r_code err)
{
	switch(err)
	{
		default:
			fprintf(stderr,"ERROR: An error occurred with return code [%hd]",err);
			break;
	}
}
