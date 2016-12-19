/*
============================================================================================
   Class for array based half-facet (AHF) data structure to store and process N-D,
   non-manifold, meshes.

   Note: Everything is indexed starting at 0!

   Also note: using a vector of structs is 2 x faster than using a vector of integers.

   Example:  2-D triangulations; in this case, a half-facet == half-edge.

   Scheme for ordering *local* topological entities:

        V2 +
           |\
           |  \
           |    \
           |      \
           |        \  E0
        E1 |          \
           |            \
           |              \
           |                \
           |                  \
        V0 +-------------------+ V1
                    E2

    EXAMPLE:  Mesh of two triangles.

        V3 +-------------------+ V2
           |\                  |
           |  \                |
           |    \        T1    |
           |      \            |
           |        \          |
           |          \        |
           |            \      |
           |     T0       \    |
           |                \  |
           |                  \|
        V0 +-------------------+ V1

   Triangle Connectivity and Sibling Half-Facet (Half-Edge) Data Struct:

   triangle |   vertices   |     sibling half-edges
    indices |  V0, V1, V2  |     E0,     E1,      E2
   ---------+--------------+-------------------------
       0    |   0,  1,  3  |  <1,1>, <NULL>,  <NULL>
       1    |   1,  2,  3  | <NULL>,  <0,0>,  <NULL>

   where <Ti,Ei> is a half-edge, where Ti is the *neighbor* triangle index, and
   Ei is the local edge index of Ti that correponds to the half-edge. <NULL> means
   there is no neighbor triangle.

   Vertex-to-Half-Edge Data Struct:

     vertex |  adjacent
    indices | half-edge
   ---------+------------
       0    |   <0,2>
       1    |   <1,2>
       2    |   <1,0>
       3    |   <0,1>

   Diagram depicting half-edges:

                   <1,0>
        V3 +-------------------+ V2
           |\                  |
           |  \          T1    |
           |    \              |
           |      \  <1,1>     |
     <0,1> |        \          | <1,2>
           |    <0,0> \        |
           |            \      |
           |              \    |
           |     T0         \  |
           |                  \|
        V0 +-------------------+ V1
                   <0,2>

   Note: in this example, only need one adjacent half-edge because there are no
         non-manifold vertices.  But we do allow for non-manifold vertices!

   Copyright (c) 12-17-2016,  Shawn W. Walker
============================================================================================
*/

#define _MESH_CC

#ifndef _BASEMESH_CC
#include "BaseMesh.cc"  // base class for all mesh topology stuff
#endif
#ifndef _BASEPTCOORD_CC
#include "BasePtCoord.cc"  // base class for all vertex coordinates
#endif

/* C++ class definition */
#define  MMC  Mesh
// template the topological and geometric dimension
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
class MMC: public BaseMesh<CELL_DIM>, public BasePtCoord<GEO_DIM>
{
public:
    MMC();
    ~MMC();

    void Clear() // clear all data
    {
        BaseMesh<CELL_DIM>::Clear();
        BasePtCoord<GEO_DIM>::Clear();
    };

private:

};

/***************************************************************************************/
/* constructor */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
MMC<CELL_DIM, GEO_DIM>::MMC() : BaseMesh<CELL_DIM>(), BasePtCoord<GEO_DIM>()
{
    // what else to do or check?
    if (GEO_DIM < CELL_DIM)
    {
        std::cout << "Desired topological dimension of a cell is " << CELL_DIM << std::endl;
        std::cout << "Desired geometric dimension is " << GEO_DIM << std::endl;
        std::cout << "Geometric dimension must be >= topological dimension!" << std::endl;
        std::exit(1);
    }
    //std::cout << "Mesh constructor..." << std::endl;
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
MMC<CELL_DIM, GEO_DIM>::~MMC()
{
    //std::cout << "Mesh destructor..." << std::endl;
}

// SWW: the only methods that go here are the ones that require BOTH cell connectivity
//      *and* vertex coordinates.

// need to "open" and "close" the cells/coordinates for appending and such...

// need adaptive refinement

// don't worry about modifying the cell data structure yet, b/c this is connected to refinement

// need topological/connectivity changes!!!!  for this, in 3-D, can think of an edge as the intersection of two half-facets?

// For my mesh generator, need to evaluate cut edges. can do this by first finding all intersected tetrahedra.  then make an edge mesh from that?
// Or just store all the unique edges of the mesh?  think of the edge mesh as a graph!  Store each edge as a pair of vertices (with smallest global index first) and the cell that the edge belongs to.

// add these methods:

/*
barycentricToCartesian - Converts the coordinates of a point from Barycentric to Cartesian
cartesianToBarycentric - Converts the coordinates of a point from Cartesian to Barycentric
barycenter             - Barycenter of triangle or tetrahedron
circumcenter           - Circumcenter of triangle or tetrahedron
faceNormal             - Triangulation face normal
featureEdges           - Triangulation sharp edges
incenter               - Incenter of triangle or tetrahedron
vertexNormal           - Triangulation vertex normal
*/

#undef MMC

/***/
