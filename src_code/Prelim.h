/*
============================================================================================
   Preliminary includes and definitions.

   Copyright (c) 12-08-2016,  Shawn W. Walker
============================================================================================
*/

#define _PRELIM_H

// std C++ classes
#include <algorithm>
#include <vector>
#include <map>
//#define NDEBUG
#include <assert.h>
#include <iostream>
#include <limits>

/* define types for indices */
typedef unsigned int SmallIndType; // type for "small" indices (i.e. 0 through ~127)
typedef unsigned int VtxIndType;   // vertex indices
typedef unsigned int CellIndType;  // cell indices, e.g. triangles, tetrahedrons

/* define NULL constants */
static const unsigned int NULL_Small = std::numeric_limits<unsigned int>::max(); // define the NULL index
static const unsigned int NULL_Vtx   = std::numeric_limits<unsigned int>::max(); // define the NULL index
static const unsigned int NULL_Cell  = std::numeric_limits<unsigned int>::max(); // define the NULL index
// we can never have an index of the maximum size (should not be a problem)

/* define data type for points */
typedef double PointType;

/***/
