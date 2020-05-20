/*
============================================================================================
   Class for array based half-facet (AHF) data structure.  This stores and processes
   multiple topological meshes (possibly non-manifold) of dimensions 1, 2, or 3.  All
   of the meshes reference a common set of vertices.  Ergo, a *single* set of point
   coordinates is also stored.  All meshes are simplex meshes.

   SWW: what about 0-dim meshes?

   Note: Everything is indexed starting at 0!

   Also note: using a vector of structs is 2 x faster than using a vector of integers.

   Copyright (c) 05-18-2020,  Shawn W. Walker
============================================================================================
*/

#define _MESH_CC

#ifndef _BASEMESH_CC
#include "BaseMesh.cc"  // base class for all mesh topology stuff
#endif
#ifndef _BASEPTCOORD_CC
#include "BasePtCoord.cc"  // base class for all vertex coordinates
#endif

/* define the mesh types we allow: None, 0-D, 1-D, 2-D, 3-D */
enum class MeshType {None, Pt, Int, Tri, Tet};

/* C++ class definition */
#define  MMC  Mesh
// template the (maximum) geometric dimension
template <SmallIndType GEO_DIM>
class MMC
{
public:
    MMC();
    MMC(SmallIndType N_0, SmallIndType N_1, SmallIndType N_2, SmallIndType N_3);
    ~MMC();
    
    /* define arrays of topologial meshes */
    std::vector< BaseMesh<0> >   PtMesh;
    std::vector< BaseMesh<1> >  IntMesh;
    std::vector< BaseMesh<2> >  TriMesh;
    std::vector< BaseMesh<3> >  TetMesh;
    /* define the set of vertex coordinates */
    BasePtCoord<GEO_DIM>  Vtx;

    // open Mesh for modification
    inline void Open()
    {
        for (std::vector< BaseMesh<0> >::iterator it=PtMesh.begin();  it!=PtMesh.end();  ++it)
            (*it).Open();
        for (std::vector< BaseMesh<1> >::iterator it=IntMesh.begin(); it!=IntMesh.end(); ++it)
            (*it).Open();
        for (std::vector< BaseMesh<2> >::iterator it=TriMesh.begin(); it!=TriMesh.end(); ++it)
            (*it).Open();
        for (std::vector< BaseMesh<3> >::iterator it=TetMesh.begin(); it!=TetMesh.end(); ++it)
            (*it).Open();
        Vtx.Open();
    };
    // close Mesh; modification is no longer allowed
    inline void Close()
    {
        for (std::vector< BaseMesh<0> >::iterator it=PtMesh.begin();  it!=PtMesh.end(); ++it)
            (*it).Close();
        for (std::vector< BaseMesh<1> >::iterator it=IntMesh.begin(); it!=IntMesh.end(); ++it)
            (*it).Close();
        for (std::vector< BaseMesh<2> >::iterator it=TriMesh.begin(); it!=TriMesh.end(); ++it)
            (*it).Close();
        for (std::vector< BaseMesh<3> >::iterator it=TetMesh.begin(); it!=TetMesh.end(); ++it)
            (*it).Close();
        Vtx.Close();
    };

    // clear all data
    void Clear()
    {
        PtMesh.clear();
        IntMesh.clear();
        TriMesh.clear();
        TetMesh.clear();
        Vtx.Clear();

        default_mesh_type  = MeshType::None; // reset
        default_mesh_index = 0;              // reset
    };
    
    // get the maximum topological dimension in the mesh
    SmallIndType Max_Topological_Dim() const
    {
        if (TetMesh.size() > 0)
            return (SmallIndType) 3;
        else if (TriMesh.size() > 0)
            return (SmallIndType) 2;
        else if (IntMesh.size() > 0)
            return (SmallIndType) 1;
        else
            return (SmallIndType) 0;
    };

private:
    // default meshes
    MeshType       default_mesh_type;
    SmallIndType   default_mesh_index;
};

/***************************************************************************************/
/* constructor */
template <SmallIndType GEO_DIM>
MMC<GEO_DIM>::MMC()
{
    MMC(0,0,0,0);
}

/***************************************************************************************/
/* constructor */
template <SmallIndType GEO_DIM>
MMC<GEO_DIM>::MMC(SmallIndType Num_0D, SmallIndType Num_1D, SmallIndType Num_2D, SmallIndType Num_3D)
{
    // start fresh
    Clear();
    
    // set default mesh dim to the maximum topological dimension
    if (Num_3D > 0)
        default_mesh_type = MeshType::Tet;
    else if (Num_2D > 0)
        default_mesh_type = MeshType::Tri;
    else if (Num_1D > 0)
        default_mesh_type = MeshType::Int;
    else if (Num_0D > 0)
        default_mesh_type = MeshType::Pt;
    else
        default_mesh_type = MeshType::None; // NULL value
    
    default_mesh_index = 0; // assume the 0th one
    
    // init the topological mesh arrays
    BaseMesh<0> Empty_PtMesh;
    PtMesh.assign(Num_0D,Empty_PtMesh);
    BaseMesh<1> Empty_IntMesh;
    IntMesh.assign(Num_1D,Empty_IntMesh);
    BaseMesh<2> Empty_TriMesh;
    TriMesh.assign(Num_2D,Empty_TriMesh);
    BaseMesh<3> Empty_TetMesh;
    TetMesh.assign(Num_3D,Empty_TetMesh);
    
    const SmallIndType MaxTopDim = Max_Topological_Dim();
    
    // what else to do or check?
    if (GEO_DIM < MaxTopDim)
    {
        std::cout << "Maximum topological dimension allowed is " << MaxTopDim << std::endl;
        std::cout << "Desired geometric dimension is " << GEO_DIM << std::endl;
        std::cout << "Geometric dimension must be >= topological dimension!" << std::endl;
        std::exit(1);
    }
    //std::cout << "Mesh constructor..." << std::endl;
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType GEO_DIM>
MMC<GEO_DIM>::~MMC()
{
    Clear();
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
barycenter             - Barycenter of triangle or tetrahedron or cell
circumcenter           - Circumcenter of triangle or tetrahedron or cell
incenter               - Incenter of triangle or tetrahedron or cell
faceNormal             - Triangulation face normal
featureEdges           - Triangulation sharp edges
vertexNormal           - Triangulation vertex normal
*/

#undef MMC

/***/
