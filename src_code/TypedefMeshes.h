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

/* define common types of meshes */
typedef Mesh<1,1> IntMesh;      // 1-D standard mesh
typedef Mesh<1,2> CurveMesh2D;  // 1-D curve mesh embedded in 2-D
typedef Mesh<1,3> CurveMesh3D;  // 1-D curve mesh embedded in 3-D
typedef Mesh<2,2> TriMesh;      // 2-D standard mesh
typedef Mesh<2,3> SurfaceMesh;  // 2-D surface mesh embedded in 3-D
typedef Mesh<3,3> TetMesh;      // 3-D standard mesh

#endif

/***/
