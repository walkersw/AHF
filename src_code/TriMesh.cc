/*
============================================================================================
   Class for array based half-facet (AHF) data structure to store and process 2-D,
   non-manifold, surface meshes.  In this case, a half-facet == half-edge.

   Note: triangles are indexed starting at 1!  Everything is indexed starting at 1!

   Also note: using a vector of structs is 2 x faster than using a vector of integers.

   Scheme for ordering *local* topological entities:

        V3 +
           |\
           |  \
           |    \
           |      \
           |        \  E1
        E2 |          \
           |            \
           |              \
           |                \
           |                  \
        V1 +-------------------+ V2
                    E3

    EXAMPLE:  Mesh of two triangles.

        V4 +-------------------+ V3
           |\                  |
           |  \                |
           |    \        T2    |
           |      \            |
           |        \          |
           |          \        |
           |            \      |
           |     T1       \    |
           |                \  |
           |                  \|
        V1 +-------------------+ V2

   Triangle Connectivity and Sibling Half-Facet (Half-Edge) Data Struct:

   triangle |   vertices   |     sibling half-edges
    indices |  V1, V2, V3  |     E1,     E2,      E3
   ---------+--------------+-------------------------
       1    |   1,  2,  4  |  <2,2>, <NULL>,  <NULL>
       2    |   2,  3,  4  | <NULL>,  <1,1>,  <NULL>

   where <ti,ei> is a half-edge, where ti is the *neighbor* triangle index, and
   ei is the local edge index of ti that correponds to the half-edge. <NULL> means
   there is no neighbor triangle.

   Vertex-to-Half-Edge Data Struct:

     vertex |  adjacent
    indices | half-edge
   ---------+------------
       1    |   <1,3>
       2    |   <2,3>
       3    |   <2,1>
       4    |   <1,2>

   Diagram depicting half-edges:

                   <2,1>
        V4 +-------------------+ V3
           |\                  |
           |  \          T2    |
           |    \              |
           |      \  <2,2>     |
     <1,2> |        \          | <2,3>
           |    <1,1> \        |
           |            \      |
           |              \    |
           |     T1         \  |
           |                  \|
        V1 +-------------------+ V2
                   <1,3>

   Note: in this example, only need one adjacent half-edge because there are no
         non-manifold vertices.  But we do allow for non-manifold vertices!

   Copyright (c) 11-05-2015,  Shawn W. Walker
============================================================================================
*/

#include "BaseMesh.cc"  // base class for all meshes

/* C++ class definition */
#define  TRI  TriMesh
class TRI: public BaseMesh<2>
{
public:
    TRI();
    ~TRI();

private:

};

/***************************************************************************************/
/* constructor */
TRI::TRI() : BaseMesh<2>()
{
    // ensure memory is clear to start
    Clear();
}

/***************************************************************************************/
/* DE-structor */
TRI::~TRI()
{
    // clear the data
    Clear();
}



// eventually, need refinement

// don't worry about modifying the structure yet, b/c this is connected to refinement

// need topological/connectivity changes!!!!  for this, in 3-D, can think of an edge as the intersection of two half-facets?

// for my mesh generator, need to evaluate cut edges. can do this by first finding all intersected tetrahedra.  then make an edge mesh from that?
// or just store all the unique edges of the mesh?  think of the edge mesh as a graph!  Store each edge as a pair of vertices (with smallest global index first) and the cell that the edge belongs to.


#undef TRI

/***/
