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
/** @file mesh.h
 *
 * @brief Definition of mesh data structures.
 * @define This library contains functions for the creation of meshes. It has been
 * designed primarily to facilitate the contruction of triangular meshes of a given
 * topological genus.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Electrical Engineering and Computer Science
 * @author Faculty of Science and Engineering
 * @author Queensland University of Technology
 * 
 * @date 10/03/2012 - 11/08/2012
 * @copyright GNU Public Lincense.
 *
 * Note: My 6-month old son was sitting on my lap while I started to write this
 *       and a few times he thought he should add a comment or two 
 *       (by bashing the keyboard) :)
 * =============================================================================
 */
 
#ifndef __MESH_H
#define __MESH_H

#include <stdio.h>
#include <string.h>
#include "vectorMath.h"
/** @brief extra memory allocted when buffers fill up.*/
#define REALLOC_EXTRA 1000
/** @brief maximum number of vertices per face.*/
#define MAX_FACE_TYPE 8
/*for ease in coordinate indexing*/
/** @brief defines X coordinate index.*/
#define X 0
/** @brief defines Y coordinate index.*/
#define Y 1
/** @brief defines Z coordinate index.*/
#define Z 2
/*return codes*/
/** @brief Return code for successful function completion.*/
#define SUCCESS 5
/** @brief Return code for successful write completion.*/
#define WRITE_SUCCESS 4
/** @brief Return code for successful read completion.*/
#define READ_SUCCESS 3
/** @brief Return code for successful deletion completion.*/
#define DELETE_SUCCESS 2
/** @brief Return code for successful insertion completion.*/
#define INSERT_SUCCESS 1
/** @brief Return code for memory failure.*/
#define OUT_OF_MEMORY 0
/** @brief Return code for failed insertion.*/
#define INSERT_FAILED -1
/** @brief Return code for failed deletion.*/
#define DELETE_FAILED -2
/** @brief Return code for failed read.*/
#define READ_FAILED -3
/** @brief Return code for failed write.*/
#define WRITE_FAILED -4
/** @brief Return code for data corruption.*/
#define DATA_CORRUPTION -126
/** @brief Return code to flag un-implemented future functionality.*/
#define NOT_IMPLEMENTED -127

/** @brief Object File Format file flag.*/
#define OFF_FORMAT 0
/** @brief Wavefront OBJ file flag.*/
#define OBJ_FORMAT 1
/** @brief STereoLithography file flag.*/
#define STL_FORMAT 2
/** @brief Virtual Reality Modeling Language file flag.*/
#define VRML_FORMAT 3



/*My Son added the following comment*/
/*pppppppppp ,m  b idc fv  df rvtngbhhhhhhyyyyyyyyyyyyyy7uuuuuuuuuuuuuuuuutyu7km*/
/** @brief A list of vertices*/ 
typedef struct vertexList_struct vertexList;
/** @brief A list faces*/ 
typedef struct faceList_struct faceList;
/** @brief A mesh*/ 
typedef struct mesh_struct mesh;
/** @brief Return code yype.*/ 
typedef short int r_code;
/** @brief The vertex list data structure.
 *  @details All vertices are packed e.g. the ith vertex is stored in 
 *  <em>verts[i*dim], \a verts[i*dim + 1],..., verts[i*dim + dim-1].</em>
 */
struct vertexList_struct 
{
    /** @brief The number of vertices in the list.*/
	int numVerts;
    /** @brief The size of the list in bytes.
     *  @note \a size = \a numVerts * \a dim .
     */
	int size;
    /** @brief The dimentionality of each vertex.*/
	int dim;
    /** @brief Memory address of array of vertex data.*/
	float *verts;/*numVerts x dim*/
};
/*Another comment by my Son*/       
/* chf */ 

/*niiiiiiiiiiiiiii uu9*/ /*again, thanks son...*/ 

/** @brief  The face list data structure. */
struct faceList_struct
{
    /** @brief The number of faces in the list.*/
	int numFaces;
    /** @brief The size of the list in bytes.
     *  @note \a size = \a numFaces * \a maxVerts .
     */
	int size;
    /** @brief The largest number of vertices a single face may have. */
	int maxVerts;
    /** @brief Memory address of array of face data.*/
	int *faces; /*numFaces x maxVertices*/
	unsigned char *faceTypes;
};

/** @brief Mesh data structure*/
struct mesh_struct
{
    /** @brief The vertex list*/
	vertexList *vList;
    /** The face list*/
	faceList *fList;
};

/*function prototypes*/

/*vertex list functions*/
vertexList * CreateVertexList(int size,int dim);
r_code InsertVertex(vertexList *vList,float *v);
r_code DeleteVertex(vertexList *vList,int i);
float * GetVertex_ptr(vertexList *vList,int i);
float * GetVertex_cp(vertexList *vList,int i);

/*face list functions*/
faceList * CreateFaceList(int size,int maxVerts);
r_code InsertFace(faceList *fList,int *face,unsigned char type);
r_code DeleteFace(faceList *fList, int i);
int * GetFace_ptr(faceList *fList,int i);
int * GetFace_cp(faceList *fList,int i);
float * FaceMidPoint(faceList *fList, int i,vertexList *vList,float *verts_out,float *m_out);
float FaceArea(faceList *fList, int i,vertexList *vList);

/*mesh functions*/
mesh * CreateMesh(int numVerts,int dim,int numFaces, int maxVerts);
mesh * CopyMesh(mesh *m);
mesh * CreateDual(mesh *m);
r_code SubDivideFaces(mesh *m);
r_code RemoveDuplicateVertices(mesh *m);
r_code CheckVertexWindings(mesh *m);
float SurfaceArea(mesh *m);
float * GeometricCentre(mesh *m);

/*Topology generation*/
mesh * CreateIcosahedron(void);
mesh * CreateTorus(void);
mesh * CreateDoubleTorus(void);
mesh * CreateMeshTopology(int numFaces,int genus);

/*I/O functions*/
mesh * LoadMesh(char *filename,unsigned char format);
r_code SaveMesh(char *filename,mesh *m, unsigned char format);

r_code WriteOFF(char *filename,mesh *m);
r_code WriteOBJ(char *filename,mesh *m);
r_code WriteSTL(char *filename,mesh *m);
r_code WriteVRML(char *filename, mesh *m);
r_code ReadOFF(char *filename,mesh **m);
r_code ReadOBJ(char *filename,mesh **m);
r_code ReadSTL(char *filename,mesh **m);
r_code ReadVRML(char *filename, mesh **m);
/*Error check codes*/
unsigned char CheckErr(r_code rc);
void PrintErrorMsg(r_code err);
#endif
