/*
============================================================================================
   Convenient typedefs for common meshes.

   Copyright (c) 12-15-2016,  Shawn W. Walker
============================================================================================
*/

#ifndef _TYPEDEFMESHES_H

#define _TYPEDEFMESHES_H

#ifndef _MESH_CC
#include "Mesh.cc"  // main general class for all meshes
#endif

/* define common types of meshes with macros */
#define IntMesh(M_Name)      Mesh<1>  M_Name(1,0,0)
#define CurveMesh2D(M_Name)  Mesh<2>  M_Name(1,0,0)
#define CurveMesh3D(M_Name)  Mesh<3>  M_Name(1,0,0)
#define TriMesh(M_Name)      Mesh<2>  M_Name(0,1,0)
#define SurfaceMesh(M_Name)  Mesh<3>  M_Name(0,1,0)
#define TetMesh(M_Name)      Mesh<3>  M_Name(0,0,1)

#endif

/***/
