/*
============================================================================================
   Convenient typedefs for common meshes.

   Copyright (c) 05-18-2020,  Shawn W. Walker
============================================================================================
*/

#ifndef _TYPEDEFMESHES_H

#define _TYPEDEFMESHES_H

#ifndef _MESH_CC
#include "Mesh.cc"  // main general class for all meshes
#endif

/* define common types of meshes with macros */
//#define PtMesh(M_Name,Vtx)       Mesh<0,1> M_Name(&(Vtx))
#define IntMesh(M_Name,Vtx)      Mesh<1,1> M_Name(&(Vtx))
#define CurveMesh2D(M_Name,Vtx)  Mesh<1,2> M_Name(&(Vtx))
#define CurveMesh3D(M_Name,Vtx)  Mesh<1,3> M_Name(&(Vtx))
#define TriMesh(M_Name,Vtx)      Mesh<2,2> M_Name(&(Vtx))
#define SurfaceMesh(M_Name,Vtx)  Mesh<2,3> M_Name(&(Vtx))
#define TetMesh(M_Name,Vtx)      Mesh<3,3> M_Name(&(Vtx))


#endif

/***/
