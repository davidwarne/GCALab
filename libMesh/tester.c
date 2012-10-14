/* File: tester.c
 *
 * Author: David J. Warne (david.warne@qut.edu.au)
 * 
 * School of Electrical Engineering and Computer Science
 * Faculty of Science and Engineering
 * Queensland University of Technology
 * 
 * Date Created: 19/03/2012
 * Last Modified: 19/03/2012
 *
 * Descritpion: test program for mesh library
 *
 * =============================================================================
 */

#include <stdio.h> 
#include <string.h>
#include "mesh.h"

int main(int argc,char **argv)
{
	r_code rc;
	int i;
	float *vert;
	int *edge;
	int *face;
	mesh *m,*m2;
	char name[255];
	if (argc == 1)
	{
		m = CreateIcosahedron();
		m2 = CreateDual(m);
		SaveMesh("icosa.off",CreateIcosahedron(),OFF_FORMAT);
		SaveMesh("torus.off",CreateTorus(),OFF_FORMAT);
		SaveMesh("icosa_180.off",CreateDual(CreateTopology(320,0)),OFF_FORMAT);
		SaveMesh("torus_180.off",CreateTopology(320,1),OFF_FORMAT);	
		SaveMesh("icosa_dual.off",m2,OFF_FORMAT);
		return 0;
	}
	rc = ReadOFF(argv[1],&m);
	for(i=0;i<(m->vList->numVerts);i++)
	{
		vert = GetVertex_ptr(m->vList,i);
		printf("%d [%f %f %f]\n",i,vert[0],vert[1],vert[2]);	
	}
/*	for(i=0;i<(m->eList->numEdges);i++)
	{
		edge = GetEdge_ptr(m->eList,i);
		printf("%d [%d %d | %d %d]\n",i,edge[V0],edge[V1],edge[F0],edge[F1]);	
	}*/
	for(i=0;i<(m->fList->numFaces);i++)
	{
		face = GetFace_ptr(m->fList,i);
		printf("%d [%d %d %d]\n",i,face[0],face[1],face[2]);	
	}
	for (i=0;i<atoi(argv[2]);i++)
		rc = SubDivideFaces(m);
	CheckErr(rc);
	sprintf(name,"%s_%d.off",argv[1],atoi(argv[2]));
	WriteOFF(name,m);
	sprintf(name,"%s_%d.obj",argv[1],atoi(argv[2]));
	WriteOBJ(name,m);
	sprintf(name,"%s_%d.stl",argv[1],atoi(argv[2]));
	WriteSTL(name,m);
	sprintf(name,"%s_%d.vrml",argv[1],atoi(argv[2]));
	WriteVRML(name,m);
	return 0;
	// here we are just testing the removal of duplicates
	m2 = CopyMesh(m);
	// we directly overwrite two verts with the second vertex
	m2->vList->verts[4*3] = m2->vList->verts[1*3]; 
	m2->vList->verts[4*3+1] = m2->vList->verts[1*3+1]; 
	m2->vList->verts[4*3+2] = m2->vList->verts[1*3+2]; 
	m2->vList->verts[8*3] = m2->vList->verts[1*3]; 
	m2->vList->verts[8*3+1] = m2->vList->verts[1*3+1]; 
	m2->vList->verts[8*3+2] = m2->vList->verts[1*3+2];

	printf("before:\n");
	for (i=0;i<m2->vList->numVerts;i++)
	{
		vert = GetVertex_ptr(m2->vList,i);
		printf("%d : %f %f %f\n",i,vert[0],vert[1],vert[2]);
	}
	for (i=0;i<m2->fList->numFaces;i++)
	{
		face = GetFace_ptr(m2->fList,i);
		printf("%d : %d %d %d\n",i,face[0],face[1],face[2]);
	}

	RemoveDuplicateVertices(m2);
	printf("after:\n");
	for (i=0;i<m2->vList->numVerts;i++)
	{
		vert = GetVertex_ptr(m2->vList,i);
		printf("%d : %f %f %f\n",i,vert[0],vert[1],vert[2]);
	}
	for (i=0;i<m2->fList->numFaces;i++)
	{
		face = GetFace_ptr(m2->fList,i);
		printf("%d : %d %d %d\n",i,face[0],face[1],face[2]);
	}
}
