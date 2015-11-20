/*
============================================================================================
   Preliminary includes and definitions.

   Copyright (c) 11-05-2015,  Shawn W. Walker
============================================================================================
*/

/* define types for indices */
typedef unsigned int SmallIndType; // type for "small" indices (i.e. 0 through ~10)
typedef unsigned int VtxIndType;   // vertex indices
typedef unsigned int CellIndType;  // cell indices, e.g. triangles, tetrahedrons

/* define constants */
static const unsigned int NULL_IND = 0; // define the NULL index

// std C++ classes
#include <algorithm>
#include <vector>
#include <map>
//#define NDEBUG
#include <assert.h>
#include <iostream>

/***/
