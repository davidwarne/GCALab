/* File: mesh.h
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 10/03/2012
 * Last Modified: 11/08/2012
 *
 * Descritpion: Definition of mesh data structures
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

#define REALLOC_EXTRA 1000
#define MAX_FACE_TYPE 8
/*for ease in coordinate indexing*/
#define X 0
#define Y 1
#define Z 2
/*return codes*/
#define SUCCESS 5
#define WRITE_SUCCESS 4
#define READ_SUCCESS 3
#define DELETE_SUCCESS 2
#define INSERT_SUCCESS 1
#define OUT_OF_MEMORY 0
#define INSERT_FAILED -1
#define DELETE_FAILED -2
#define READ_FAILED -3
#define WRITE_FAILED -4
#define DATA_CORRUPTION -126
#define NOT_IMPLEMENTED -127

#define OFF_FORMAT 0
#define OBJ_FORMAT 1
#define STL_FORMAT 2
#define VRML_FORMAT 3



/*My Son added the following comment*/
/*pppppppppp ,m  b idc fv  df rvtngbhhhhhhyyyyyyyyyyyyyy7uuuuuuuuuuuuuuuuutyu7km*/
typedef struct vertexList_struct vertexList;
typedef struct faceList_struct faceList;
typedef struct mesh_struct mesh;
typedef short int r_code;
/* vertex list data structure, all vertices are packed 
 *  e.g. the ith vertex is stored in verts[i*dim],verts[i*dim + 1],
 *  ...,verts[i*dim + dim-1]
 */
struct vertexList_struct 
{
	int numVerts;
	int size;
	int dim;
	float *verts;/*numVerts x dim*/
};
/*Another comment by my Son*/       
/* chf */ 

/*niiiiiiiiiiiiiii uu9*/ /*again, thanks son...*/ 

/* face list data structure */
struct faceList_struct
{
	int numFaces;
	int size;
	int maxVerts;
	int *faces; /*numFaces x maxVertices*/
	unsigned char *faceTypes;
};

/*mesh data structure*/
struct mesh_struct
{
	vertexList *vList;
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
float * FaceMidPoint(faceList *fList, int i,vertexList *vList);
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
mesh * CreateTopology(int numFaces,int genus);

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
