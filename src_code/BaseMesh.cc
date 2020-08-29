/*
============================================================================================
   Base class for array based half-facet (AHF) data structure to store and process meshes.

   Note: this base class is used in deriving the 0-D, 1-D, 2-D, and 3-D mesh classes,
   as well as arbitrarily higher dimensions (all simplex meshes).
   Note: no vertex coordinates are stored in this class; this is purely topological.
   Note: everything is indexed from 0!

   Scheme for ordering *local* topological entities
   ------------------------------------------------

   RULE: local "facet" Fi is always *opposite* local vertex Vi.

   DIM=1. Interval: facets are vertices (points)

        V0 +-------------------+ V1

   i.e. "facet" F0 = V1 (which is opposite V0), and "facet" F1 = V0 (which is opposite V1)

   DIM=2. Triangle: facets are edges (line segments)

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

   i.e. "facet" F0 = E0 (which is opposite V0), "facet" F1 = E1 (which is opposite V1), etc...

   DIM=3: Tetrahedron: facets are faces (triangles)

                 V3 +
                   /|\
                  / | \
                 |  |  \
                |   |   \
               |    |    \
               |    |  F1 \               F0 (opposite V0)
              | F2  |      \
              |     |       \
             |   V0 +--------+ V2
             |     /      __/
            |    /  F3 __/
            |  /    __/
           | /   __/
          |/  __/
       V1 +--/

   i.e. facet F0 = [V1, V2, V3] (which is opposite V0),
        facet F1 = [V0, V3, V2] (which is opposite V1),
        facet F2 = [V0, V1, V3] (which is opposite V2),
        facet F3 = [V0, V2, V1] (which is opposite V3).

   Higher DIM: the pattern continues...

   EXAMPLE:  Mesh of two triangles.  In this case, a half-facet == half-edge.

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

   Copyright (c) 05-18-2020,  Shawn W. Walker
============================================================================================
*/

#define _BASEMESH_CC

#ifndef _PRELIM_H
#include "Prelim.h" // basic typedefs and includes
#endif
#ifndef _VTX2HALFFACET_MAPPING_CC
#include "Vtx2HalfFacet_Mapping.cc" // class for handling the Vertex-to-Half-Facet Mapping
#endif

#ifndef _BASICSTRUCTS_H
#include "BasicStructs.h" // basic structs for the BaseMesh class
#endif

/* C++ class definition */
#define  BM  BaseMesh
// template the cell topological dimension
// Note: the number of cell facets equals top. dim. + 1
template <SmallIndType CELL_DIM>
class BM
{
public:
    BM();
    ~BM();
    // open the mesh for modification
    inline void Open() { Mesh_Open = true; };
    // close the mesh; modification is no longer allowed
    inline void Close() { Mesh_Open = false; };
    // generic check for if mesh is open for modification
    inline bool Is_Mesh_Open() const
    {
        if (!Mesh_Open)
        {
            std::cout << "Mesh is not open for modification!" << std::endl;
            std::cout << "     You must first use the 'Open' method." << std::endl;
        }
        return Mesh_Open;
    };

    void Clear() // clear all data
    {
        Cell.clear();
        Vtx2HalfFacets.Clear();
        v2hfs.Clear();

        Mesh_Open = true; // re-open the mesh for modification
    };

    // get the topological dimension
    inline SmallIndType Top_Dim() const { return CELL_DIM; };
    // allocate room for specified number of elements
    void Reserve_Cells(const CellIndType&);
    // set the cell data with a given array
    void Set_Cell_Data(const VtxIndType*, const CellIndType&);
    // append a bunch of cells (in a given array) to the end of the array
    void Append_Cell_Data(const VtxIndType*, const CellIndType&);
    // append one cell to the end of the array
    void Append_Cell(const VtxIndType*);
    void Append_Cell(const VtxIndType&, const VtxIndType&);
    void Append_Cell(const VtxIndType&, const VtxIndType&, const VtxIndType&);
    void Append_Cell(const VtxIndType&, const VtxIndType&, const VtxIndType&, const VtxIndType&);
    void Append_Cell_0D(const VtxIndType&);
    // ... and update the intermediate data structure v2hfs (build incrementally)
    void Append_Cell_And_Update(const VtxIndType*);
    void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&);
    void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&, const VtxIndType&);
    void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&, const VtxIndType&, const VtxIndType&);
    void Append_Cell_0D_And_Update(const VtxIndType&);
    // get number of elements stored
    inline CellIndType Num_Cells() const { return (CellIndType) Cell.size(); };

    // retrieve one cell (writeable)
    inline VtxIndType* Get_Cell_vtx(const CellIndType&);
    // this one is for the neighbor query
    inline HalfFacetType* Get_Cell_halffacet(const CellIndType&);
    inline void Get_Cell(const CellIndType&, VtxIndType*, HalfFacetType*);
    // retrieve one cell in its struct form (writeable)
    //inline CellType& Get_Cell(const CellIndType&);
    inline CellSimplexType<CELL_DIM>& Get_Cell(const CellIndType&);

    // retrieve one cell (read-only)
    inline const VtxIndType* Get_Cell_vtx(const CellIndType&) const;
    // this one is for the neighbor query
    inline const HalfFacetType* Get_Cell_halffacet(const CellIndType&) const;
    inline void Get_Cell(const CellIndType&, const VtxIndType*&, const HalfFacetType*&) const;
    // retrieve one cell in its struct form (read-only)
    //inline const CellType& Get_Cell(const CellIndType&) const;
    inline const CellSimplexType<CELL_DIM>& Get_Cell(const CellIndType&) const;

    // get unique set of vertices
    void Get_Unique_Vertices(std::vector<VtxIndType>&);
    void Display_Unique_Vertices();
    // get number of vertices referenced in Cell
    VtxIndType Num_Vtx()
    {
        std::vector<VtxIndType> uv;
        Get_Unique_Vertices(uv);
        return (VtxIndType) uv.size();
    };
    // get maximum vertex index referenced in Cell
    VtxIndType Max_Vtx_Index()
    {
        std::vector<VtxIndType> uv;
        Get_Unique_Vertices(uv);
        return (VtxIndType) uv.back();
    };

    // map vertex indices to new indices
    void Reindex_Vertices(const VtxIndType&, const VtxIndType*);
    
    // finalize the data structures for determining mesh connectivity
    //    (i.e. neighbors, vtx2half-facet mapping, etc.) and *close* the mesh.
    void Finalize_Mesh_Connectivity();

    /* all public routines below this need the mesh to be finalized to output
           correct information, i.e. the mesh should be "closed" and all internal
           data structures updated. This is done by building the sibling half-facet
           structure, and filling out the Vtx2HalfFacets mapping. All of this is
           automatically done by the "Finalize_Mesh_Connectivity" method. */

    // return const reference to v2hfs (an intermediate, internal data structure)
    const Vtx2HalfFacet_Mapping& Get_v2hfs() const { const Vtx2HalfFacet_Mapping& c_v2hfs = v2hfs; return c_v2hfs; };
    // get read access to Vtx2HalfFacets
    inline const Vtx2HalfFacet_Mapping& Get_Vtx2HalfFacets() const { return Vtx2HalfFacets; };

    //void Display_Unique_Vertices() const { Vtx2HalfFacets.Display_Unique_Vertices(); };

    // returns a unique set of all edges of the mesh
    void Get_Edges(std::vector<MeshEdgeType>&) const;
    // print out all the edges of the mesh
    void Display_Edges() const;
    // test if a pair of vertices is connected by an edge
    bool Is_Connected(const VtxIndType&, const VtxIndType&) const;
    bool Is_Connected(const VtxIndType vv[2]) const;
    // returns all cell indices attached to a given edge
    void Get_Cells_Attached_To_Edge(const MeshEdgeType&, std::vector<CellIndType>&) const;
    // returns all half-facets that are referenced by only one cell;
    //         i.e. the half-facets that are on the boundary of the mesh
    void Get_FreeBoundary(std::vector<HalfFacetType>&) const;

    // print out cell connectivity and sibling half-facet data
    void Display_Cell(const CellIndType& ci=NULL_Cell) const;
    // print out half-facets attached to vertex in intermediate data structure
    void Display_v2hfs(const VtxIndType& vi=NULL_Vtx) const;
    // print out half-facets attached to vertex for final data structure
    void Display_Vtx2HalfFacets(const VtxIndType& vi=NULL_Vtx) const;
    // returns all cell indices attached to a given vertex
    void Get_Cells_Attached_To_Vertex(const VtxIndType&, std::vector<CellIndType>&) const;
    // returns all cell indices attached to a given vertex and facet-connected to the given cell
    void Get_Cells_Attached_To_Vertex(const VtxIndType&, const CellIndType&, std::vector<CellIndType>&) const;
    // display cell indices attached to given vertex
    void Display_Cells_Attached_To_Vertex(const VtxIndType&) const;
    // determine if two cells (that share a vertex) are connected by a sequence of facet neighbors
    bool Two_Cells_Are_Facet_Connected(const VtxIndType&, const CellIndType&, const CellIndType&) const;
    // display if two cells (that both contain the given vertex) are facet-sequence-connected
    void Display_Two_Cells_Are_Facet_Connected(const VtxIndType&, const CellIndType&, const CellIndType&) const;
    // get all half-facets attached to a given half-facet
    void Get_HalfFacets_Attached_To_HalfFacet(const HalfFacetType&, std::vector<HalfFacetType>&) const;
    // display attached half-facets
    void Display_HalfFacets_Attached_To_HalfFacet(const HalfFacetType&) const;
    // get a unique set of non-manifold half-facets
    void Get_Nonmanifold_HalfFacets(std::vector<HalfFacetType>&) const;
    // display all non-manifold half-facets
    void Display_Nonmanifold_HalfFacets() const;
    // get the set of non-manifold vertices
    void Get_Nonmanifold_Vertices(std::vector<VtxIndType>&) const;
    // display all non-manifold vertices
    void Display_Nonmanifold_Vertices() const;

    /* these methods here are for debugging purposes only.
       the casual user should never use them! */
    // note: these call internal private versions.
    //       see the private methods (with similar names) for more info.
    void Finalize_v2hfs_DEBUG(bool bfs=true);
    void Build_Sibling_HalfFacets_DEBUG();
    void Build_Vtx2HalfFacets_DEBUG();

protected:
    /* main data storage */
    typedef CellSimplexType<CELL_DIM> CellSimplex_DIM; // convenient
    // connectivity and sibling half-facet data
    std::vector<CellSimplex_DIM>  Cell;
    // referenced vertices in Cell and (possibly multiple) attached half-facet(s)
    Vtx2HalfFacet_Mapping         Vtx2HalfFacets;

    double       Cell_Reserve_Buffer; // amount of extra memory to allocate when re-allocating
                                      // Cell and Vtx2HalfFacets (number between 0.0 and 1.0).
    VtxIndType   Estimate_Size_Vtx2HalfFacets; // estimate of the size to allocate in Vtx2HalfFacets.

    // intermediate data structure for building sibling half-facet information
    Vtx2HalfFacet_Mapping   v2hfs; // Note: for a given vertex, it references multiple half-facets.
    // Note: this data structure will NOT NECESSARILY store all referenced vertices
    //       in the triangulation.  This is because the vertex with smallest index
    //       will never be referenced (for example).  This is an internal structure that
    //       is only used to construct the sibling half-facet information (stored in Cell).

    /* flag to indicate if mesh cells may be added or modified.
       true  = cells can be added, modified
       false = the mesh cells cannot be changed! */
    bool Mesh_Open;

    // internal cell struct access
    inline CellSimplex_DIM& Get_Cell_struct(const CellIndType&); // for writeable
    inline const CellSimplex_DIM& Get_Cell_struct(const CellIndType&) const; // for read-only

    /* private methods for building the mesh connectivity,
           i.e. sibling half-facets, vtx2half-facet mapping, etc.) */
    // finalize intermediate data structure
    // (input = true means build from scratch using cell connectivity data)
    void Finalize_v2hfs(bool bfs=true);
    // finishing filling out the sibling half-facet data structure
    void Build_Sibling_HalfFacets();
    // build the final Vtx2HalfFacets data struct
    void Build_Vtx2HalfFacets();

    // append half facets of given cell to v2hfs
    void Append_Half_Facets(const CellIndType&, const VtxIndType v[(CELL_DIM+1)]);
    // mapping used for generating the sibling half-facets
    void Vtx2Adjacent(const VtxIndType&, const CellIndType&, const SmallIndType&, VtxIndType* va) const;

    // amount extra to reserve when finding a variable number of cells attached to a vertex
    static const SmallIndType cell_attach_chunk = 5*CELL_DIM;
    // recursion routines for public methods above
    void Get_Cells_Attached_To_Vertex_Recurse(const VtxIndType&, const CellIndType&, std::vector<CellIndType>&) const;
    // recursive call for public method above
    bool Two_Cells_Are_Facet_Connected_Recurse(const VtxIndType&, const CellIndType&, const CellIndType&,
                                               const CellIndType&, CellIndType&) const;

    /* baby routines */
    // given a global vertex, this determines which local vertex it is in the cell
    inline SmallIndType Get_Local_Vertex_Index_In_Cell(const VtxIndType&, const CellSimplex_DIM&) const;
    // map from local vertex in cell to local facets in cell that contain that vertex
    inline void Get_Local_Facets_Sharing_Local_Vertex(const SmallIndType& vi, SmallIndType*) const;
    // map local facet index to local vertices contained in that facet
    inline void Get_Local_Vertices_Of_Local_Facet(const SmallIndType&, SmallIndType*) const;

    // get the global vertices in a cell's facet
    inline void Get_Global_Vertices_In_Facet
          (const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi, VtxIndType* fv) const;
    // get the vertices in a facet that are adjacent to the given vertex
    inline void Get_Adj_Vertices_In_Facet
          (const VtxIndType*, const VtxIndType&, VtxIndType*) const;
    // simple max operation
    VtxIndType Get_Vertex_With_Largest_Index_In_Facet(const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi) const;
    // determine if adjacent vertices are equal
    inline bool Adj_Vertices_In_Facet_Equal(const VtxIndType* a, const VtxIndType* b) const;

private:
    // internal helper routine
    void Push_Back_Cell(const CellSimplex_DIM&);
};

/***************************************************************************************/
/* constructor */
template <SmallIndType CELL_DIM>
BM<CELL_DIM>::BM ()
{
    Mesh_Open = true; // the mesh starts out as open for modification

    Cell_Reserve_Buffer = 0.2; // allocate an extra 20% when re-allocating
    Estimate_Size_Vtx2HalfFacets = 0;
    // ensure memory is clear to start
    Clear();

    //std::cout << "BaseMesh constructor..." << std::endl;
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType CELL_DIM>
BM<CELL_DIM>::~BM()
{
    // clear the data
    Clear();

    //std::cout << "BaseMesh destructor..." << std::endl;
}

/***************************************************************************************/
/* Allocate memory to hold triangulation of given size (plus a little). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Reserve_Cells(const CellIndType& Num_C)
{
    if (!Is_Mesh_Open())
        return;

    // compute the actual size to allocate for the cells
    const CellIndType Actual_Cell_SIZE = (CellIndType) ((1.0 + Cell_Reserve_Buffer) * Num_C);
    Cell.reserve(Actual_Cell_SIZE);
    // guess on what to reserve for the intermediate data structure
    v2hfs.Reserve((VtxIndType) (CELL_DIM+1) * Num_C);
}

/***************************************************************************************/
/* Set the cell data all at once.
   Note: input is the global vertex indices of the corners of all the cell simplices. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Set_Cell_Data(const VtxIndType* Cell_Data, const CellIndType& Num_Cell_Data)
{
    // allocate
    Reserve_Cells(Num_Cell_Data);
    // init to NULL
    CellSimplex_DIM CL;
    CL.Clear();
    Cell.assign(Num_Cell_Data,CL);
    
    // now fill it in
    for (CellIndType ii = 0; ii < Num_Cell_Data; ++ii)
    {
        Cell[ii].Set(Cell_Data + (CELL_DIM+1)*ii);
    }
}

/***************************************************************************************/
/* Internal helper routine for pushing a cell to the end of the list. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Push_Back_Cell(const CellSimplex_DIM& CL)
{
    if (!Is_Mesh_Open())
        return;

    // first check that there is room; if not, then reserve more space
    const CellIndType Size_Est = (CellIndType) Cell.size();
    if (Size_Est >= (CellIndType) Cell.capacity())
        Reserve_Cells(Size_Est+2); // then out of room, so reserve more space
    Cell.push_back(CL);
}

/***************************************************************************************/
/* Append a cell element to the end of the list.
   Note: input is the global vertex indices of the corners of all the cell simplices. */
// append a bunch of cells at once
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_Data(const VtxIndType* Cell_Data, const CellIndType& Num_Cell_Data)
{
    const CellIndType New_Total_Cells = Num_Cells() + Num_Cell_Data;
    Reserve_Cells(New_Total_Cells);
    
    // now append them
    for (CellIndType ii = 0; ii < Num_Cell_Data; ++ii)
    {
        CellSimplex_DIM  CL;
        CL.Set(Cell_Data + (CELL_DIM+1)*ii);
        Push_Back_Cell(CL);
    }
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType* vtx)
{
    CellSimplex_DIM CL;
    CL.Set(vtx);
    /* append cell to the end */
    Push_Back_Cell(CL);
}
// 1-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType& v0, const VtxIndType& v1)
{
    CellSimplex_DIM CL;
    CL.Set(v0,v1);
    /* append cell to the end */
    Push_Back_Cell(CL);
}
// 2-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType& v0, const VtxIndType& v1, const VtxIndType& v2)
{
    CellSimplex_DIM CL;
    CL.Set(v0,v1,v2);
    /* append cell to the end */
    Push_Back_Cell(CL);
}
// 3-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType& v0, const VtxIndType& v1,
                               const VtxIndType& v2, const VtxIndType& v3)
{
    CellSimplex_DIM CL;
    CL.Set(v0,v1,v2,v3);
    /* append cell to the end */
    Push_Back_Cell(CL);
}
// 0-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_0D(const VtxIndType& v0)
{
    CellSimplex_DIM CL;
    CL.Set_0D(v0);
    /* append cell to the end */
    Push_Back_Cell(CL);
}

/***************************************************************************************/
/* Append a cell element to the end of the list,
   and build the intermediate v2hfs structure (incrementally). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_And_Update(const VtxIndType* vtx)
{
    // get the next cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Cell(vtx);

    // now "ci" is the *current* cell index
    Append_Half_Facets(ci, vtx);
}
// 1-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_And_Update(const VtxIndType& v0, const VtxIndType& v1)
{
    // get the next cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Cell(v0,v1);

    // now "ci" is the *current* cell index
    const VtxIndType vtx[2] = {v0, v1};
    Append_Half_Facets(ci, vtx);
}
// 2-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_And_Update(const VtxIndType& v0, const VtxIndType& v1, const VtxIndType& v2)
{
    // get the next cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Cell(v0,v1,v2);

    // now "ci" is the *current* cell index
    const VtxIndType vtx[3] = {v0, v1, v2};
    Append_Half_Facets(ci, vtx);
}
// 3-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_And_Update(const VtxIndType& v0, const VtxIndType& v1,
                                          const VtxIndType& v2, const VtxIndType& v3)
{
    // get the next cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Cell(v0,v1,v2,v3);

    // now "ci" is the *current* cell index
    const VtxIndType vtx[4] = {v0, v1, v2, v3};
    Append_Half_Facets(ci, vtx);
}
// 0-D
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_0D_And_Update(const VtxIndType& v0)
{
    // get the next cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Cell_0D(v0);

    // now "ci" is the *current* cell index
    const VtxIndType vtx[1] = {v0};
    Append_Half_Facets(ci, vtx);
}

/***************************************************************************************/
/* Get parts (or all) of a specific cell's data (given the cell index). */
// WRITEABLE
template <SmallIndType CELL_DIM>
inline VtxIndType* BM<CELL_DIM>::Get_Cell_vtx(const CellIndType& ci) // writeable
{
    if (!Is_Mesh_Open())
    {
        std::cerr << "Fatal error in 'Get_Cell_vtx'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci].vtx;
}
template <SmallIndType CELL_DIM>
inline HalfFacetType* BM<CELL_DIM>::Get_Cell_halffacet(const CellIndType& ci) // writeable
{
    if (!Is_Mesh_Open())
    {
        std::cerr << "Fatal error in 'Get_Cell_halffacet'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci].halffacet;
}
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Cell(const CellIndType& ci, VtxIndType* vtx_p, HalfFacetType* hf_p) // writeable
{
    if (!Is_Mesh_Open())
    {
        std::cerr << "Fatal error in 'Get_Cell'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    vtx_p = Cell[ci].vtx;
    hf_p  = Cell[ci].halffacet;
}
/* Get pointer to specific mesh cell data (given the cell index). */
template <SmallIndType CELL_DIM>
//inline CellType& BM<CELL_DIM>::Get_Cell(const CellIndType& ci) // writeable
inline CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell(const CellIndType& ci) // writeable
{
    if (!Is_Mesh_Open())
    {
        std::cerr << "Fatal error in 'Get_Cell'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci];
}
// READ-only
template <SmallIndType CELL_DIM>
inline const VtxIndType* BM<CELL_DIM>::Get_Cell_vtx(const CellIndType& ci) const // read-only
{
    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci].vtx;
}
template <SmallIndType CELL_DIM>
inline const HalfFacetType* BM<CELL_DIM>::Get_Cell_halffacet(const CellIndType& ci) const // read-only
{
    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci].halffacet;
}
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Cell(const CellIndType& ci, const VtxIndType*& vtx_p, const HalfFacetType*& hf_p) const // read-only
{
    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    vtx_p = Cell[ci].vtx;
    hf_p  = Cell[ci].halffacet;
}
/* Get pointer to specific mesh cell data (given the cell index). */
template <SmallIndType CELL_DIM>
//inline const CellType& BM<CELL_DIM>::Get_Cell(const CellIndType& ci) const // read-only
inline const CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell(const CellIndType& ci) const // read-only
{
    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci];
}

/***************************************************************************************/
/* get unique list of vertices (the non-const version). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Unique_Vertices(std::vector<VtxIndType>& uv)
{
    uv.clear();
    
    if (Mesh_Open)
    {
        const CellIndType NC = Cell.size(); 
        const SmallIndType Num_Local_Vtx = (SmallIndType) (CELL_DIM+1);
        // make it big enough to include disjoint cells
        uv.reserve( Num_Local_Vtx * NC + 1);

        // add all vertices of each cell
        for (CellIndType ci = 0; ci < NC; ++ci)
        {
            const VtxIndType* c_vtx = Cell[ci].vtx;
            // loop through each vertex of the current (simplex) cell
            for (SmallIndType vi = 0; vi < CELL_DIM+1; ++vi)
            {
                const VtxIndType global_vi = c_vtx[vi];
                uv.push_back(global_vi);
            }
        }
        std::sort(uv.begin(),uv.end());
        
        // get unique set of vertices
        std::vector<VtxIndType>::iterator it_end;
        it_end = std::unique_copy(uv.begin(), uv.end(), uv.begin());
        // get number of unique vertices
        const VtxIndType LENGTH = (unsigned int) std::distance(uv.begin(),it_end);
        // resize
        uv.resize(LENGTH);
    }
    else
    {
        Vtx2HalfFacets.Get_Unique_Vertices(uv);
    }
}

/***************************************************************************************/
/* print unique list of vertices */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Unique_Vertices()
{
    if (Mesh_Open)
    {
        // extract all the vertex indices
        std::vector<VtxIndType> unique_vertices;
        Get_Unique_Vertices(unique_vertices);

        // now print them out
        std::cout << "Unique list of vertex indices: " << std::endl;
        std::cout << *(unique_vertices.begin());
        for (std::vector<VtxIndType>::const_iterator vi=unique_vertices.begin()+1; vi!=unique_vertices.end(); ++vi)
        {
            std::cout << ", " << (*vi);
        }
        std::cout << std::endl;
    }
    else
    {
        Vtx2HalfFacets.Display_Unique_Vertices();
    }
}

/***************************************************************************************/
/* re-index the vertices in the mesh.
   example:  new_index = new_indices[old_index] */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Reindex_Vertices(const VtxIndType& num_new_indices, const VtxIndType* new_indices)
{
    // the mesh must be *open* to do this.
    if (!Is_Mesh_Open())
        return;

    // basic check
    if (num_new_indices < Max_Vtx_Index())
    {
        std::cerr << "Fatal error in 'BaseMesh.Reindex_Vertices'!" << std::endl;
        std::cerr << "    The given list of indices is shorter than the max vertex index" << std::endl;
        std::cerr << "    referenced by cells in the mesh." << std::endl;
        std::exit(1);
    }
    
    // go through all the cells, and map those vertices
    const CellIndType NC = Num_Cells();
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        VtxIndType* c_vtx = Cell[ci].vtx;
        // loop through each vertex of the current (simplex) cell
        for (SmallIndType vi = 0; vi < CELL_DIM+1; ++vi)
        {
            c_vtx[vi] = new_indices[c_vtx[vi]];
        }
    }
    
    // go through Vtx2HalfFacets and map its vertices
    std::vector<VtxHalfFacetType>& V2HF = Vtx2HalfFacets.Get_VtxMap();
    for (std::vector<VtxHalfFacetType>::iterator it=V2HF.begin(); it!=V2HF.end(); ++it)
        {
            (*it).vtx = new_indices[(*it).vtx];
        }
}

/***************************************************************************************/
/* finalize the data structures for determining mesh connectivity, i.e.
   determine neighbors (sibling half-facets), vtx2half-facet mapping, etc.)
   and *close* the mesh. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Finalize_Mesh_Connectivity()
{
    // the mesh must be *open* to do this.
    if (!Is_Mesh_Open())
        return;

    // this sequence of commands must be used!
    Finalize_v2hfs(true); // setup an intermediate structure
    Build_Sibling_HalfFacets();
    Build_Vtx2HalfFacets();

    // now *close* the mesh to further modification
    Close();
}

/***************************************************************************************/
/* fill "edges" with all the edges in the mesh
   Note: this is a unique list of *sorted* edges, i.e. each edge [v0, v1] satisfies
         v0 < v1, where v0, v1 are the global vertex indices of the edge end points.
         Moreover, the edges are in ascending order. */
bool Mesh_Edge_Ascending_order(const MeshEdgeType& Ea, const MeshEdgeType& Eb)
{
    if (Ea.vtx[0]==Eb.vtx[0])
        return (Ea.vtx[1] < Eb.vtx[1]);
    else
        return (Ea.vtx[0] < Eb.vtx[0]);
}
bool Mesh_Edge_equality(const MeshEdgeType& Ea, const MeshEdgeType& Eb)
{
    return Ea.Equal(Eb);
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Edges(std::vector<MeshEdgeType>& edges) const
{
    const CellIndType NC = Num_Cells();
    const SmallIndType Num_Local_Edge = (SmallIndType) ((CELL_DIM+1) * CELL_DIM / 2);
    // make it big enough to include repeats
    edges.reserve( Num_Local_Edge * NC + 1);

    // add all edges of each cell
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        // loop through each local edge of the current (simplex) cell
        // i.e. loop through all the distinct pairs of vertices in the cell
        for (SmallIndType vi = 0; vi < CELL_DIM+1; ++vi)
        {
            MeshEdgeType EE;
            const VtxIndType* c_vtx = Get_Cell_vtx(ci);
            for (SmallIndType vj = vi+1; vj < CELL_DIM+1; ++vj)
            {
                EE.Set_Sorted(c_vtx[vi], c_vtx[vj]);
                edges.push_back(EE);
            }
        }
    }

    // sort all the edges
    std::sort(edges.begin(), edges.end(), Mesh_Edge_Ascending_order);
    // now get a unique list (no repeats!)
    std::vector<MeshEdgeType>::iterator it;
    it = std::unique(edges.begin(), edges.end(), Mesh_Edge_equality);
    edges.resize(std::distance(edges.begin(),it));
}

/***************************************************************************************/
/* print out all the edges in the mesh (this uses "Get_Edges"). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Edges() const
{
    std::vector<MeshEdgeType> edges;
    Get_Edges(edges);
    if (edges.size() > 0)
    {
        // print all mesh edges
        std::cout << "Display all edges of the mesh (a unique list):" << std::endl;
        std::cout << "Edge:  [ vertex #0, vertex #1 ]" << std::endl;
        for (std::vector<MeshEdgeType>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
        {
            (*it).Print();
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout << "There are NO mesh edges!" << std::endl;
    }
    std::cout << std::endl;
}

/***************************************************************************************/
/* test if a pair of vertices is connected by an edge. */
template <SmallIndType CELL_DIM>
bool BM<CELL_DIM>::Is_Connected(const VtxIndType& v0, const VtxIndType& v1) const
{
    if (v0==v1)
    {
        std::cout << "Input vertices v0, v1 are the *same*!" << std::endl;
        return true; // trivial
    }

    // get all cells attached to v0
    std::vector<CellIndType> attached_cells;
    Get_Cells_Attached_To_Vertex(v0, attached_cells);
    // check if v1 is in any of the cells
    for (std::vector<CellIndType>::const_iterator it=attached_cells.begin(); it!=attached_cells.end(); ++it)
    {
        const VtxIndType* cell_vtx;
        cell_vtx = Get_Cell_vtx(*it);
        for (SmallIndType ii = 0; ii < CELL_DIM+1; ++ii)
        {
            if (cell_vtx[ii]==v1)
                return true; // found it!
        }
    }

    return false; // was not found
}
template <SmallIndType CELL_DIM>
bool BM<CELL_DIM>::Is_Connected(const VtxIndType vv[2]) const
{
    return Is_Connected(vv[0], vv[1]);
}

/***************************************************************************************/
/* returns all cell indices attached to a given edge. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Cells_Attached_To_Edge(const MeshEdgeType& EE, std::vector<CellIndType>& attached_cells) const
{
    const VtxIndType& v0 = EE.vtx[0];
    const VtxIndType& v1 = EE.vtx[1];

    // get all cells attached to v0
    std::vector<CellIndType> attached_to_v0;
    Get_Cells_Attached_To_Vertex(v0, attached_to_v0);
    std::sort(attached_to_v0.begin(),attached_to_v0.end());

    // get all cells attached to v1
    std::vector<CellIndType> attached_to_v1;
    Get_Cells_Attached_To_Vertex(v1, attached_to_v1);
    std::sort(attached_to_v1.begin(),attached_to_v1.end());

    // take the intersection
    attached_cells.resize(attached_to_v0.size() + attached_to_v1.size(), NULL_Cell); // init size
    std::vector<CellIndType>::iterator it;
    it = std::set_intersection(attached_to_v0.begin(), attached_to_v0.end(),
                               attached_to_v1.begin(), attached_to_v1.end(), attached_cells.begin());
    attached_cells.resize(it - attached_cells.begin());
}

/***************************************************************************************/
/* returns all half-facets that are referenced by only one cell;
   i.e. the half-facets that are on the boundary of the mesh. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_FreeBoundary(std::vector<HalfFacetType>& bdy) const
{
    const CellIndType NC = Num_Cells();
    const SmallIndType Num_Local_Facets = CELL_DIM+1;
    // reserve some space to start
    bdy.reserve( (CellIndType) (Num_Local_Facets * NC / 4) + 10 );

    // check all cells
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        const HalfFacetType* Cell_HFs = Get_Cell_halffacet(ci);
        // loop through each local facet of the current (simplex) cell
        for (SmallIndType fi = 0; fi < CELL_DIM+1; ++fi)
        {
            HalfFacetType HF;
            if (Cell_HFs[fi].Is_Null())
            {
                // this facet has no neighbor!
                HF.Set(ci, fi); // so store this bdy facet

                const CellIndType current_capacity = bdy.capacity();
                // make sure there is room to store
                if (bdy.size() >= current_capacity)
                    bdy.reserve( (CellIndType) 1.5*current_capacity );
                bdy.push_back(HF);
            }
        }
    }
}

/***************************************************************************************/
/* print cell connectivity and sibling haf-facets. "ci" is the index of a
   specific cell; if ci=NULL_Cell, then print all cells. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Cell(const CellIndType& main_ci) const
{
    if (main_ci==NULL_Cell)
    {
        // then print all the cells
        std::cout << "Display connectivity of all cells:" << std::endl;
        std::cout << "Cell #    |        Vertices        |   Sibling Half-Facets" << std::endl;
        for (CellIndType ci = 0; ci < Num_Cells(); ++ci)
        {
            const CellSimplex_DIM& CL = Get_Cell_struct(ci);
            std::cout << ci << "  |  " << CL.vtx[0]; // print cell # and first vertex
            for (SmallIndType kk = 1; kk < (CELL_DIM+1); ++kk)
                std::cout << ", " << CL.vtx[kk];
            std::cout << "   |   ";
            CL.halffacet[0].Print();
            for (SmallIndType kk = 1; kk < (CELL_DIM+1); ++kk)
            {
                std::cout << ", ";
                CL.halffacet[kk].Print();
            }
            std::cout << std::endl;
        }
    }
    else
    {
        // then print ONE cell
        std::cout << "Display connectivity of cell #" << main_ci << ":" << std::endl;
        std::cout << "        Vertices        |   Sibling Half-Facets" << std::endl;

        const CellSimplex_DIM& CL = Get_Cell_struct(main_ci);
        std::cout << "  " << CL.vtx[0]; // print first vertex
        for (SmallIndType kk = 1; kk < (CELL_DIM+1); ++kk)
            std::cout << ", " << CL.vtx[kk];
        std::cout << "   |   ";
        CL.halffacet[0].Print();
        for (SmallIndType kk = 1; kk < (CELL_DIM+1); ++kk)
        {
            std::cout << ", ";
            CL.halffacet[kk].Print();
        }
        std::cout << std::endl;
    }
}

/***************************************************************************************/
/* print (multiple) half-facets attached to a given vertex
   (from intermediate data structure). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_v2hfs(const VtxIndType& vi) const
{
    std::cout << "'v2hfs':" << std::endl;
    v2hfs.Display_Half_Facets(vi);
}

/***************************************************************************************/
/* print half-facets attached to a given vertex
   (for final data structure). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Vtx2HalfFacets(const VtxIndType& vi) const
{
    std::cout << "'Vtx2HalfFacets':" << std::endl;
    Vtx2HalfFacets.Display_Half_Facets(vi);
}

/***************************************************************************************/
/* returns all cell indices (in "cell_array") that are attached to vertex "vi".
   Note: this requires the sibling half-facet data (in the Cell variable) to be
   built before this can be used; also, Vtx2HalfFacets must be finished as well. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Cells_Attached_To_Vertex(const VtxIndType& vi, std::vector<CellIndType>& cell_array) const
{
    // init the array
    cell_array.clear();
    if (vi==NULL_Vtx) return; // the vertex is null, so do nothing.

    // get the attached half-facets
    std::pair <std::vector<VtxHalfFacetType>::const_iterator,
               std::vector<VtxHalfFacetType>::const_iterator> RR;
    const MedIndType Num_HF = Vtx2HalfFacets.Get_Half_Facets(vi, RR);

    // guess the number of cells attached to vi
    cell_array.reserve(Num_HF*cell_attach_chunk);

    // loop through each half-facet
    // (each one corresponds to a connected component of the mesh)
    for (std::vector<VtxHalfFacetType>::const_iterator it=RR.first; it!=RR.second; ++it)
    {
        std::vector<CellIndType> temp_array;
        Get_Cells_Attached_To_Vertex(vi, (*it).ci, temp_array);
        // store the found cells in cell_array
        const CellIndType old_size = (CellIndType) cell_array.size();
        cell_array.resize(old_size + temp_array.size(), NULL_Cell);
        std::copy(temp_array.begin(), temp_array.end(), cell_array.begin()+old_size );
    }
}

/***************************************************************************************/
/* returns all cell indices (in "cell_array") that are attached to vertex "vi" and are
   facet-connected to the cell "ci". Note: this requires the sibling half-facet data
   (in the Cell variable) to be built before this can be used. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Cells_Attached_To_Vertex(const VtxIndType& vi, const CellIndType& ci, std::vector<CellIndType>& cell_array) const
{
    // init the array
    cell_array.clear();
    if ( (vi==NULL_Vtx) || (ci==NULL_Cell) ) return; // the vertex or starting cell is null, so do nothing.

    // guess the number of cells attached to vi
    cell_array.reserve(cell_attach_chunk);

    // do a recursive search to collect all the cells
    Get_Cells_Attached_To_Vertex_Recurse(vi, ci, cell_array);
}
/* recursive function call for the above method. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Cells_Attached_To_Vertex_Recurse(const VtxIndType& vi, const CellIndType& ci, std::vector<CellIndType>& cell_array) const
{
    // if we have been searching too long...
    if (cell_array.size() > 100000)
    {
        // then quit!
        std::cout << "ERROR in 'Get_Cells_Attached_To_Vertex'..." << std::endl;
        std::cout << "    Recursion depth is too great." << std::endl;
        std::cout << "    There should not be more than 100,000 cells attached to a single vertex!" << std::endl;
        return;
    }

    // cell is null, so nothing left to do
    if (ci==NULL_Cell) return; // this can happen if the neighbor cell does not exist

    std::vector<CellIndType>::iterator it;
    it = std::find(cell_array.begin(), cell_array.end(), ci); // see if ci is already in the array
    // if ci is already in the array
    if (it!=cell_array.end())
        return; // then we have already visited this cell, so done
    else
    {
        // access the cell
        const CellSimplex_DIM& CL = Get_Cell_struct(ci);
        // check again that the vertex is actually in the cell
        const SmallIndType local_vi = Get_Local_Vertex_Index_In_Cell(vi, CL);
        if (local_vi==NULL_Small) return; // cell does not actually contain the vertex (this should not happen)

        // add the cell to the list
        const CellIndType current_size = (CellIndType) cell_array.size();
        if (current_size == cell_array.capacity())
            cell_array.reserve(current_size + cell_attach_chunk);
        cell_array.push_back(ci);

        if (CELL_DIM > 0)
        {
            // get the local facets that share the vertex
            SmallIndType* local_facet = new SmallIndType[CELL_DIM];
            Get_Local_Facets_Sharing_Local_Vertex(local_vi, local_facet);

            // loop through the facets
            for (SmallIndType fi = 0; fi < CELL_DIM; ++fi)
                if (local_facet[fi]!=NULL_Small)
                {
                    // determine the neighbor cell on the "other side" of that facet and search it
                    Get_Cells_Attached_To_Vertex_Recurse(vi, CL.halffacet[local_facet[fi]].ci, cell_array);
                }
                // else the neighbor does not exist, so do nothing
            delete(local_facet);
        }
    }
}

/***************************************************************************************/
/* print all cell indices that are attached to vertex "vi".
   This prints out the separate connected components.
   Note: this requires the sibling half-facet data (in the Cell variable) to be
   built before this can be used; also, Vtx2HalfFacets must be finished as well. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Cells_Attached_To_Vertex(const VtxIndType& vi) const
{
    if (vi==NULL_Vtx)
    {
        std::cout << "Vertex is invalid, so display nothing!" << std::endl;
        return; // the vertex is null, so do nothing.
    }

    // get the attached half-facets
    std::pair <std::vector<VtxHalfFacetType>::const_iterator,
               std::vector<VtxHalfFacetType>::const_iterator> RR;
    Vtx2HalfFacets.Get_Half_Facets(vi, RR);

    std::cout << "Vertex #" << vi << " is attached to the following cells:" << std::endl;

    // loop through each half-facet
    // (each one corresponds to a connected component of the mesh)
    SmallIndType COUNT = 0;
    for (std::vector<VtxHalfFacetType>::const_iterator it=RR.first; it!=RR.second; ++it)
    {
        ++COUNT;
        std::cout << "component #" << COUNT << ": ";
        std::vector<CellIndType> cell_array;
        Get_Cells_Attached_To_Vertex(vi, (*it).ci, cell_array);
        std::cout << cell_array[0];
        // print it
        for (CellIndType cc=1; cc < cell_array.size(); ++cc)
            std::cout << ", " << cell_array[cc];
        std::cout << std::endl;
    }
}

/***************************************************************************************/
/* this returns true/false if two cells (that share the same vertex) are connected
   by a "chain" of half-facet neighbors. Note: this is useful for determining when
   two cells are in the same "connected component" of the mesh (this is important
   when the mesh is *not* a manifold). */
template <SmallIndType CELL_DIM>
bool BM<CELL_DIM>::Two_Cells_Are_Facet_Connected(const VtxIndType& vi, const CellIndType& ci_a, const CellIndType& ci_b) const
{
    if ( (vi==NULL_Vtx) || (ci_a==NULL_Cell) || (ci_b==NULL_Cell) ) return false; // something is null, so do nothing.

    // init recursion depth count
    CellIndType Depth_Count = 0;
    // init current cell and target cell
    const CellIndType   Start_Cell = ci_a;
    const CellIndType Current_Cell = ci_a;
    const CellIndType  Target_Cell = ci_b;

    // do a recursive search to find the target
    return Two_Cells_Are_Facet_Connected_Recurse(vi, Start_Cell, Current_Cell, Target_Cell, Depth_Count);
}
/* recursive function call for the above method. */
template <SmallIndType CELL_DIM>
bool BM<CELL_DIM>::Two_Cells_Are_Facet_Connected_Recurse(const VtxIndType& vi, const CellIndType& start, const CellIndType& current,
                                                         const CellIndType& target, CellIndType& Depth_Count) const
{
    // if we have been searching too long...
    if (Depth_Count > 100000)
    {
        // then quit!
        std::cout << "ERROR in 'Two_Cells_Are_Facet_Connected'..." << std::endl;
        std::cout << "    Recursion depth is too great." << std::endl;
        std::cout << "    There should not be more than 100,000 cells attached to a single vertex!" << std::endl;
        return false;
    }

    // if loop back to the beginning...
    if (Depth_Count > 0)
        if (current==start) return false; // we did not find the target

    // current cell is null, so nothing left to do
    if (current==NULL_Cell) return false; // this can happen if the neighbor cell does not exist

    // if current matches the target
    if (current==target)
        return true; // then we found it!
    else
    {
        // access the cell
        const CellSimplex_DIM& CL = Get_Cell_struct(current);
        // check again that the vertex is actually in the cell
        const SmallIndType local_vi = Get_Local_Vertex_Index_In_Cell(vi, CL);
        if (local_vi==NULL_Small) return false; // cell does not actually contain the vertex (this should not happen)

        // get the local facets that share the vertex
        SmallIndType* local_facet = new SmallIndType[CELL_DIM];
        Get_Local_Facets_Sharing_Local_Vertex(local_vi, local_facet);

        // keep counting
        Depth_Count++;

        // loop through the facets
        bool CONNECTED = false; // init
        for (SmallIndType fi = 0; fi < CELL_DIM; ++fi)
            if (local_facet[fi]!=NULL_Small)
            {
                // determine the neighbor cell on the "other side" of that facet and search it
                CONNECTED = Two_Cells_Are_Facet_Connected_Recurse(
                                           vi, start, CL.halffacet[local_facet[fi]].ci, target, Depth_Count);
                if (CONNECTED)
                    break; // it was found!
            }
            // else the neighbor does not exist, so do nothing

        delete(local_facet);
        return CONNECTED;
    }
}

/***************************************************************************************/
/* print information stating whether two cells (that share the same vertex) are
   facet-connected. Note: this is useful for determining when
   two cells are in the same "connected component" of the mesh (this is important
   when the mesh is *not* a manifold).
*/
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Two_Cells_Are_Facet_Connected
                  (const VtxIndType& vi, const CellIndType& ci_a, const CellIndType& ci_b) const
{
    if (vi==NULL_Vtx)
    {
        std::cout << "Vertex is invalid, so display nothing!" << std::endl;
        return; // the vertex is null, so do nothing.
    }
    if ( (ci_a==NULL_Cell) || (ci_b==NULL_Cell) )
    {
        std::cout << "One of the cells is invalid, so display nothing!" << std::endl;
        return; // cell is null, so do nothing.
    }

    const bool CONNECTED = Two_Cells_Are_Facet_Connected(vi, ci_a, ci_b);
    //cout << "Check if two cells are facet connected AND share a vertex:" << endl;
    if (CONNECTED)
    {
        std::cout << "Cell #" << ci_a << " and Cell #" << ci_b << " are facet-connected *and* share vertex #" << vi << "." << std::endl;
    }
    else
    {
        // make sure both cells contain the vertex
        bool contain_a = false;
        const CellSimplex_DIM CL_a = Get_Cell_struct(ci_a);
        for (SmallIndType kk = 0; kk < (CELL_DIM+1); ++kk)
            if (CL_a.vtx[kk]==vi)
            {
                contain_a = true;
                break;
            }
        bool contain_b = false;
        const CellSimplex_DIM CL_b = Get_Cell_struct(ci_b);
        for (SmallIndType kk = 0; kk < (CELL_DIM+1); ++kk)
            if (CL_b.vtx[kk]==vi)
            {
                contain_b = true;
                break;
            }
        if ( contain_a && contain_b )
        {
            std::cout << "Cell #" << ci_a << " and Cell #" << ci_b << " both share vertex #" << vi << " but are *not* facet connected." << std::endl;
        }
        else
        {
            std::cout << "Cell #" << ci_a << " and Cell #" << ci_b << " do *not* both share vertex #" << vi << "." << std::endl;
        }
    }
}

/***************************************************************************************/
/* Get all half-facets attached to a given half-facet.  Note that all of these attached
   half-facets refer to the *same* geometrically defined facet in the mesh.

   The output of this method is a std::vector.
   Note: the output also contains the given half-facet.
   Note: this routine requires the sibling half-facet data to be built first. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_HalfFacets_Attached_To_HalfFacet
                  (const HalfFacetType& hf_in, std::vector<HalfFacetType>& attached_hf) const
{
    attached_hf.clear(); // start fresh
    attached_hf.reserve(2); // most of the time, there are only 2 cells sharing a facet,
                            // unless the mesh is not a manifold!

    // verify that the given half-facet is not NULL
    if (hf_in.Is_Null()) return;

    // put in the initial half-facet
    attached_hf.push_back(hf_in);

    // cycle through all the neighbors and store them
    MedIndType COUNT = 0;
    while (COUNT < 100000) // allow for up to 100,000 neighbors!
    {
        ++COUNT;
        // get the next half-facet
        const HalfFacetType& current_hf = attached_hf.back();
        const CellSimplex_DIM CL = Get_Cell_struct(current_hf.ci);
        const HalfFacetType& next_hf = CL.halffacet[current_hf.fi];

        // if the neighbor does not exist, then stop!
        if (next_hf.Is_Null()) break;

        // if we get back to the starting half-facet, stop!
        if (next_hf.Equal(hf_in)) break;
        // else store it
        const MedIndType current_capacity = attached_hf.capacity();
        // reserve in chunks
        if (attached_hf.size()==current_capacity)
            attached_hf.reserve(current_capacity + 5);
        attached_hf.push_back(next_hf);
    }

    if (COUNT >= 100000)
    {
        // then quit!
        std::cout << "ERROR in 'Get_HalfFacets_Attached_To_HalfFacet'..." << std::endl;
        std::cout << "    Number of neighbors is too large." << std::endl;
        std::cout << "    There should not be more than 100,000 cells attached to a single facet!" << std::endl;
    }
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_HalfFacets_Attached_To_HalfFacet(const HalfFacetType& hf_in) const
{
    std::vector<HalfFacetType> attached;
    Get_HalfFacets_Attached_To_HalfFacet(hf_in, attached);

    std::cout << "The half-facets attached to <" << hf_in.ci << "," << hf_in.fi << "> are:" << std::endl;
    for (MedIndType jj = 0; jj < attached.size(); ++jj)
        std::cout << "<" << attached[jj].ci << "," << attached[jj].fi << ">" << std::endl;
}

/***************************************************************************************/
/* get a unique set of non-manifold half-facets. Note: this requires the sibling
   half-facet data to be completed.  This returns a vector of half-facets, each
   defining a *distinct* non-manifold half-facet. */
bool non_manifold_sort_function(HalfFacetType A, HalfFacetType B) { return (A.ci < B.ci); };
bool non_manifold_equal_function(HalfFacetType& A, HalfFacetType& B) { return ( A.Equal(B) ); };
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Nonmanifold_HalfFacets(std::vector<HalfFacetType>& non_manifold_hf) const
{
    non_manifold_hf.clear(); // start fresh
    non_manifold_hf.reserve(10); // init

    // go thru all the elements
    const CellIndType NC = Num_Cells();
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        const CellSimplex_DIM CL = Get_Cell_struct(ci);
        // loop through all the (half) facets
        for (SmallIndType fi = 0; fi < (CELL_DIM+1); ++fi)
        {
            // get the facet neighbor
            const HalfFacetType& neighbor_hf = CL.halffacet[fi];
            const SmallIndType& n_ci = neighbor_hf.ci;
            const SmallIndType& n_fi = neighbor_hf.fi;

            // if the neighbor is not-NULL
            if (!neighbor_hf.Is_Null())
            {
                const CellSimplex_DIM N_CL = Get_Cell_struct(n_ci);
                // if the neighbor half-facet looks back at the cell we started at
                if (N_CL.halffacet[n_fi].ci==ci)
                {
                    // and if the local facet does *not* match where we started
                    if (N_CL.halffacet[n_fi].fi!=fi)
                    {
                        // then: two distinct facets of a cell are joined together!!
                        // this should not happen!
                        std::cout << "ERROR in 'Get_Nonmanifold_HalfFacets':" << std::endl;
                        std::cout << "      Two facets of a cell are siblings; this should not happen!" << std::endl;
                        assert(N_CL.halffacet[n_fi].fi==fi);
                    }
                    // else the local facet matches where we started,
                    //      so the starting half-facet only has one neighbor,
                    //      i.e. it is *manifold* (do nothing).
                }
                else // there is more than one neighbor, so this half-facet is *not* manifold
                {
                    // get all of the neighbors (including the starting half-facet)
                    std::vector<HalfFacetType> vec_neighbor;
                    Get_HalfFacets_Attached_To_HalfFacet(neighbor_hf, vec_neighbor);

                    // get the neighbor half-facet with the largest cell index
                    const HalfFacetType& MAX_hf = *std::max_element(vec_neighbor.begin(), vec_neighbor.end(), non_manifold_sort_function);

                    // store that half-facet
                    const MedIndType current_capacity = non_manifold_hf.capacity();
                    // make sure there is room to store
                    if (non_manifold_hf.size()==current_capacity)
                        non_manifold_hf.reserve(2*current_capacity);
                    non_manifold_hf.push_back(MAX_hf);
                }
            }
            // else there is no neighbor so this half-facet is *manifold* (do nothing)
        }
    }

    // now clean it up by removing duplicate half-facets
    std::sort(non_manifold_hf.begin(), non_manifold_hf.end(), non_manifold_sort_function);
    std::vector<HalfFacetType>::iterator it;
    it = std::unique(non_manifold_hf.begin(), non_manifold_hf.end(), non_manifold_equal_function);
    non_manifold_hf.resize( (MedIndType) std::distance(non_manifold_hf.begin(), it) );
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Nonmanifold_HalfFacets() const
{
    std::vector<HalfFacetType> non_manifold_hf;
    Get_Nonmanifold_HalfFacets(non_manifold_hf);

    const MedIndType NUM = non_manifold_hf.size();
    if (NUM==0)
        std::cout << "There are *no* non-manifold half-facets." << std::endl;
    else // there is at least 1
    {
        std::cout << "These are all the non-manifold half-facets in the mesh (output: <cell index, local facet index>):" << std::endl;
        for (MedIndType jj = 0; jj < NUM; ++jj)
            std::cout << "<" << non_manifold_hf[jj].ci << "," << non_manifold_hf[jj].fi << ">" << std::endl;
    }
}

/***************************************************************************************/
/* get the set of non-manifold vertices. Note: this requires the 'Vtx2HalfFacets'
   variable to be filled in.  This returns a vector of vertex indices, each
   defining a *distinct* non-manifold vertex. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Nonmanifold_Vertices(std::vector<VtxIndType>& non_manifold_vtx) const
{
    non_manifold_vtx.clear(); // start fresh
    non_manifold_vtx.reserve(10); // init

    // access the vtx-to-half-facet data
    const std::vector<VtxHalfFacetType>& V2HF = Vtx2HalfFacets.Get_VtxMap();

    // check everything (note: V2HF is already sorted.)
    // Note: no need to check the last entry.
    for(std::vector<VtxHalfFacetType>::const_iterator it = V2HF.begin(); it!=V2HF.end()-1; ++it)
    {
        std::vector<VtxHalfFacetType>::const_iterator next_it = it+1;
        const VtxIndType& current_vtx = (*it).vtx;
        const VtxIndType&    next_vtx = (*next_it).vtx;
        // if a vertex shows up more than once, then it is a *non-manifold* vertex
        if (current_vtx==next_vtx)
        {
            // add this vertex to the output vector
            const VtxIndType current_capacity = non_manifold_vtx.capacity();
            // make sure there is room to store
            if (non_manifold_vtx.size()==current_capacity)
                non_manifold_vtx.reserve(2*current_capacity);

            non_manifold_vtx.push_back(current_vtx);
        }
    }

    // now clean it up by removing duplicate vertices
    std::sort(non_manifold_vtx.begin(), non_manifold_vtx.end());
    std::vector<VtxIndType>::iterator stop;
    stop = std::unique(non_manifold_vtx.begin(), non_manifold_vtx.end());
    non_manifold_vtx.resize( (VtxIndType) std::distance(non_manifold_vtx.begin(), stop) );
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Nonmanifold_Vertices() const
{
    std::vector<VtxIndType> non_manifold_vtx;
    Get_Nonmanifold_Vertices(non_manifold_vtx);

    const VtxIndType NUM = non_manifold_vtx.size();
    if (NUM==0)
        std::cout << "There are *no* non-manifold vertices." << std::endl;
    else // there is more than 1
    {
        std::cout << "These are all the non-manifold vertex indices:" << std::endl;
        for (VtxIndType jj = 0; jj < NUM; ++jj)
            std::cout << non_manifold_vtx[jj] << std::endl;
    }
}

/***************************************************************************************/
/* see 'Finalize_v2hfs' for more info. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Finalize_v2hfs_DEBUG(bool bfs)
{
    Finalize_v2hfs(bfs);
}
/***************************************************************************************/
/* see 'Build_Sibling_HalfFacets' for more info. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Build_Sibling_HalfFacets_DEBUG()
{
    Build_Sibling_HalfFacets();
}
/***************************************************************************************/
/* see 'Build_Vtx2HalfFacets' for more info. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Build_Vtx2HalfFacets_DEBUG()
{
    Build_Vtx2HalfFacets();
}

/***************************************************************************************/
/* Get reference to specific mesh cell data (given the cell index). */
template <SmallIndType CELL_DIM>
inline CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell_struct(const CellIndType& ci)
{
    if (!Is_Mesh_Open())
    {
        std::cerr << "Fatal error in 'Get_Cell_struct'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci];
}
template <SmallIndType CELL_DIM>
inline const CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell_struct(const CellIndType& ci) const // read-only
{
    // ci must be in [0, Num_Cells), and not invalid
    assert((ci < Num_Cells()) && (ci!=NULL_Cell));

    return Cell[ci];
}

/***************************************************************************************/
/* store adjacent half-facets to vertices in intermediate data structure. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Finalize_v2hfs(bool Build_From_Scratch)
{
    if (!Is_Mesh_Open())
        return;

    if (Build_From_Scratch)
    {
        v2hfs.Clear(); // start fresh
        // record all vertex indices (with duplicates)
        const CellIndType NC = Num_Cells();
        v2hfs.Reserve((VtxIndType) (CELL_DIM+1) * NC);
        for (CellIndType ci = 0; ci < NC; ++ci)
        {
            const CellSimplex_DIM& CL = Cell[ci];
            Append_Half_Facets(ci, CL.vtx);
        }
    }
    // don't forget to sort!
    v2hfs.Sort();
}

/***************************************************************************************/
/* fill in the sibling half-facet data structure.
   Note: this updates the internal data of "Cell". */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Build_Sibling_HalfFacets()
{
    if ( (!Is_Mesh_Open()) || (CELL_DIM==0) )
        return;

    // need some temporary variables
    SmallIndType arr_size;
    if (CELL_DIM > 0)
        arr_size = CELL_DIM-1;
    else
        arr_size = 0;
    VtxIndType* Adj_Vtx       = new VtxIndType[arr_size];
    VtxIndType* Adj_Vtx_First = new VtxIndType[arr_size];
    VtxIndType* Adj_Vtx_Other = new VtxIndType[arr_size];

    // go thru all the elements
    const CellIndType NC = Num_Cells();
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        // loop through all the facets
        for (SmallIndType ff = 0; ff < (CELL_DIM+1); ++ff)
        {
            if (Cell[ci].halffacet[ff].Is_Null()) // i.e. it is uninitialized
            {
                // find the vertex with largest ID in the face ff
                const VtxIndType MaxVtx = Get_Vertex_With_Largest_Index_In_Facet(Cell[ci].vtx, ff);
                // get vertices in the facet ff of the current element that is adjacent to MaxVtx
                Vtx2Adjacent(MaxVtx, ci, ff, Adj_Vtx); // see below...

                // find all half-facets that are attached to MaxVtx
                std::pair<std::vector<VtxHalfFacetType>::const_iterator,
                          std::vector<VtxHalfFacetType>::const_iterator>  RR; // make range variable
                const MedIndType Num_HF = v2hfs.Get_Half_Facets(MaxVtx, RR);

                // update sibling half-facets in Cell to be a cyclic mapping...
                if (Num_HF > 0) // then there is at least one half-facet attached to MaxVtx
                {
                    // keep track of consecutive pairs in cyclic mapping
                    std::vector<VtxHalfFacetType>::const_iterator  Current, Next;

                    // Note: all the attached half-facets are attached to MaxVtx.
                    //       we say a half-facet is *valid* if its adjacent vertices match the
                    //       adjacent vertices in the facet of the original cell (... see above).

                    // find the first valid half-facet...
                    std::vector<VtxHalfFacetType>::const_iterator  Start = RR.second-1; // default value
                    for (std::vector<VtxHalfFacetType>::const_iterator it=RR.first; it!=RR.second; ++it)
                    {
                        // get vertices in the half-facet that is adjacent to MaxVtx
                        Vtx2Adjacent(MaxVtx, (*it).ci, (*it).fi, Adj_Vtx_First);

                        // if the adjacent facet vertices match, then this half-facet is valid
                        if (Adj_Vertices_In_Facet_Equal(Adj_Vtx_First, Adj_Vtx))
                        {
                            // ... and save it
                            Start = it;
                            break;
                        }
                    }
                    // init Current to Start
                    Current = Start;

                    // loop through the remaining half-facets
                    for (Next=Current+1; Next!=RR.second; ++Next)
                    {
                        // get vertices in the half-facet that is adjacent to MaxVtx
                        Vtx2Adjacent(MaxVtx, (*Next).ci, (*Next).fi, Adj_Vtx_Other);

                        // if the half-facet is valid
                        if (Adj_Vertices_In_Facet_Equal(Adj_Vtx_Other, Adj_Vtx))
                        {
                            // in the Current cell and half-facet, write the Next half-facet
                            CellSimplex_DIM& Current_Cell = Get_Cell_struct((*Current).ci);
                            Current_Cell.halffacet[(*Current).fi].ci = (*Next).ci;
                            Current_Cell.halffacet[(*Current).fi].fi = (*Next).fi;
                            // Update Current to Next
                            Current = Next;
                        }
                    }
                    // don't forget to close the cycle:
                    // if Current is different from Start
                    if (Current!=Start) // i.e. it cannot refer to itself!
                    {
                        // then in the Current cell and half-facet, write the Start half-facet
                        CellSimplex_DIM& Current_Cell = Get_Cell_struct((*Current).ci);
                        Current_Cell.halffacet[(*Current).fi].ci = (*Start).ci;
                        Current_Cell.halffacet[(*Current).fi].fi = (*Start).fi;
                    }
                    // Note: Current and Start are guaranteed to be valid at this point!
                }
                else
                {
                    // error check!
                    std::cout << "ERROR: nothing is attached to the largest index vertex in the facet!" << std::endl;
                    HalfFacetType  hf;
                    v2hfs.Get_Half_Facet(MaxVtx, hf);
                    assert(!hf.Is_Null()); // this should stop the program
                }
            }
        }
    }
    // clear temp variables
    delete(Adj_Vtx_Other);
    delete(Adj_Vtx_First);
    delete(Adj_Vtx);

    // now that we no longer need v2hfs, we can delete it
    Estimate_Size_Vtx2HalfFacets = v2hfs.Size(); // for estimating size of Vtx2HalfFacets
    v2hfs.Clear();
}

/***************************************************************************************/
/* Build the final vertex-to-adjacent half-facets data structure.
   In general, a vertex is attached to (or contained in) many half-facets.  But we
   only need to store one of the half-facets (for each vertex), because the rest can
   be found by a simple local search using the Sibling Half-Facets.  If the vertex
   is a NON-manifold vertex, then we need to store more than one half-facet in order
   to be able to easily find *all* half-facets (and all cells) that contain that
   vertex.  The number of half-facets needed (for a single vertex) depends on
   how "non-manifold" the vertex is. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Build_Vtx2HalfFacets()
{
    if (!Is_Mesh_Open())
        return;

    Vtx2HalfFacets.Clear(); // start fresh
    // make guess on how much space to allocate
    Vtx2HalfFacets.Reserve(Estimate_Size_Vtx2HalfFacets);

    // allocate a temp variable and initialize to all false,
    //    because we have not visited any half-facets yet.
    typedef CellSimplexMarkedType<CELL_DIM> CellMarked_DIM;
    CellMarked_DIM mm;
    for (SmallIndType kk = 0; kk < (CELL_DIM+1); ++kk)
        mm.facet[kk] = false;
    const CellIndType NC = Num_Cells();
    std::vector<CellMarked_DIM> Marked(NC,mm);

    // store one-to-one mapping from local vertices to local facets
    SmallIndType lv_to_lf[CELL_DIM+1];
    if (CELL_DIM==1)
    {
        lv_to_lf[0] = 0;
        lv_to_lf[1] = 1;
    }
    else
    {
        for (SmallIndType kk = 0; kk < CELL_DIM; ++kk)
            lv_to_lf[kk] = kk+1;
        lv_to_lf[CELL_DIM] = 0;
    }

    // loop through each cell
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        // loop through each local vertex of the current cell
        for (SmallIndType local_vi = 0; local_vi < (CELL_DIM+1); ++local_vi)
        {
            // current half-facet is <ci,local_fi>, local_fi = lv_to_lf(local_vi)
            const SmallIndType local_fi = lv_to_lf[local_vi];
            if (!Marked[ci].facet[local_fi])
            {
                // store this half-facet
                VtxHalfFacetType vhf;
                const VtxIndType Global_Vtx = Cell[ci].vtx[local_vi]; // get global vertex
                vhf.vtx = Global_Vtx;
                vhf.ci  = ci;
                vhf.fi  = local_fi;
                Vtx2HalfFacets.Append(vhf);
                // denote that we have visited it!
                Marked[ci].facet[local_fi] = true;

                // get the cells attached to the current vertex, that are also
                //     facet-connected to the current cell
                std::vector<CellIndType> attached_cells;
                Get_Cells_Attached_To_Vertex(Global_Vtx, ci, attached_cells);
                // loop thru the attached cells and mark corresponding half-facets as also visited
                for (std::vector<CellIndType>::iterator it = attached_cells.begin(); it < attached_cells.end(); ++it)
                {
                    const CellIndType ci_hat = *it;
                    // the local vertex index within ci_hat
                    const SmallIndType local_vi_hat = Get_Local_Vertex_Index_In_Cell(Global_Vtx, Cell[ci_hat]);
                    // corresponding local facet index
                    const SmallIndType local_fi_hat = lv_to_lf[local_vi_hat];
                    // mark it!
                    Marked[ci_hat].facet[local_fi_hat] = true;
                }
            }
        }
    }

    Vtx2HalfFacets.Sort(); // now this data structure is usable!

    /* Give border half-facets (i.e. half-facets with no siblings) higher priority.
       This allows us to easily identify vertices that are on the boundary of the mesh. */

    // loop through each cell
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        // loop through each local facet of the current cell
        for (SmallIndType local_fi = 0; local_fi < (CELL_DIM+1); ++local_fi)
        {
            // if this half-facet has no sibling
            if (Cell[ci].halffacet[local_fi].Is_Null())
            {
                if (CELL_DIM > 0)
                {
                    // get the local vertices of the local facet
                    SmallIndType* local_vtx = new SmallIndType[CELL_DIM];
                    Get_Local_Vertices_Of_Local_Facet(local_fi, local_vtx);

                    // for each vertex of the current half-facet
                    for (SmallIndType local_vi = 0; local_vi < CELL_DIM; ++local_vi)
                    {
                        // get the global vertex
                        const VtxIndType global_vi = Cell[ci].vtx[ local_vtx[local_vi] ];

                        // get the half-facets attached to the vertex
                        std::pair <std::vector<VtxHalfFacetType>::iterator,
                                   std::vector<VtxHalfFacetType>::iterator> RR;
                        const MedIndType Num_HF = Vtx2HalfFacets.Get_Half_Facets(global_vi, RR);

                        // find the half-facet in the connected component corresponding to ci,
                        //      and replace it with the half-facet with no sibling.
                        if (Num_HF==0)
                        {
                            // error!
                            std::cout << "ERROR: the first part of 'Build_Vtx2HalfFacets' missed this vertex: "
                                      << global_vi << "." << std::endl;
                            HalfFacetType temp_hf;
                            Vtx2HalfFacets.Get_Half_Facet(global_vi, temp_hf);
                            assert(!temp_hf.Is_Null()); // this should stop the program
                        }
                        else if (Num_HF==1)
                        {
                            // in this case, it is obvious what to replace
                            (*(RR.first)).ci = ci;
                            (*(RR.first)).fi = local_fi;
                        }
                        else
                        {
                            // there is more than one connected component,
                            //    so we need to find the correct one to replace.
                            for (std::vector<VtxHalfFacetType>::iterator hf_it = RR.first; hf_it != RR.second; ++hf_it)
                            {
                                const bool CONNECTED = Two_Cells_Are_Facet_Connected(global_vi, ci, (*hf_it).ci);
                                if (CONNECTED)
                                {
                                    // then replace this one
                                    (*hf_it).ci = ci;
                                    (*hf_it).fi = local_fi;
                                    break; // can stop looking
                                }
                            }
                        }
                    }
                    delete(local_vtx);
                }
            }
        }
    }
    // Note: we don't have to sort again,
    //       because the half-facets are ordered by the attached vertex
}

/***************************************************************************************/
/* Append half-facets to v2hfs struct.
   ci = cell index, vtx = array of vertex indices of the cell, ci. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Half_Facets(const CellIndType& ci, const VtxIndType vtx[(CELL_DIM+1)])
{
    VtxHalfFacetType vhf;
    vhf.ci = ci; // store the cell index

    /* loop through all the half-facets of the cell */
    for (SmallIndType fi = 0; fi < (CELL_DIM+1); ++fi)
    {
        // associate (local #fi) half-facet with the vertex with largest index
        //           within that half-facet
        vhf.fi = fi;
        const VtxIndType VV = Get_Vertex_With_Largest_Index_In_Facet(vtx, fi);
        v2hfs.Append(VV, vhf);
    }
}

/***************************************************************************************/
/* Given a vertex and half-facet, return the other vertices in the half-facet.
   If the half-facet does not contain the given vertex, then the output is NULL.
   Note: <ci,fi> is a half-facet, where ci is a cell index, and fi is a local
   facet index. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Vtx2Adjacent(const VtxIndType& v_in, const CellIndType& ci, const SmallIndType& fi,
                                VtxIndType* v_adj) const
{
    if ( (fi > CELL_DIM) ) std::cout << "ERROR: local facet index is invalid!" << std::endl;
    assert( (fi <= CELL_DIM) );

    if (CELL_DIM >= 2)
    {
        // get the vertices in the half-facet "fi"
        VtxIndType* vtx_in_hf = new VtxIndType[CELL_DIM];
        const CellSimplex_DIM CL = Get_Cell_struct(ci);
        Get_Global_Vertices_In_Facet(CL.vtx, fi, vtx_in_hf);

        // now return the vertices in the half-facet, EXCEPT v_in (i.e. the adjacent vertices)
        Get_Adj_Vertices_In_Facet(vtx_in_hf, v_in, v_adj);
        delete(vtx_in_hf);
    }
    // else do nothing!  There are no adjacent vertices.
}

/***************************************************************************************/
/* Find the *local* index of the given *global* vertex within the given cell. */
template <SmallIndType CELL_DIM>
inline SmallIndType BM<CELL_DIM>::Get_Local_Vertex_Index_In_Cell(const VtxIndType& vi, const CellSimplex_DIM& CL) const
{
    for (SmallIndType ii=0; ii < (CELL_DIM+1); ++ii)
        if (CL.vtx[ii]==vi) return ii;

    return NULL_Small;
}

/***************************************************************************************/
/* return the local facet indices that share a given local vertex.
   Note: indexing starts at 0. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Local_Facets_Sharing_Local_Vertex(const SmallIndType& vi, SmallIndType* facet) const
{
    if ( (vi > CELL_DIM) ) std::cout << "ERROR: vertex index vi is too small or too large!" << std::endl;
    assert( (vi <= CELL_DIM) );

    // note: vertex vi is opposite facet vi
    //       so, facet vi does NOT contain vertex vi

    // get the local facets that contain local vertex vi
    for (SmallIndType kk = 0; kk < vi; ++kk)
        facet[kk] = kk;
    for (SmallIndType kk = vi+1; kk < (CELL_DIM+1); ++kk)
        facet[kk-1] = kk;
}

/***************************************************************************************/
/* Return the local vertices that are attached to a given local facet.
   Note: indexing starts at 0. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Local_Vertices_Of_Local_Facet(const SmallIndType& fi, SmallIndType* vert) const
{
    // hahah, we can reuse this!
    Get_Local_Facets_Sharing_Local_Vertex(fi, vert);
}

/***************************************************************************************/
/* given a cell's vertices and a local facet, return the global vertices contained
   in that facet. Note: fi indexes from 0. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Global_Vertices_In_Facet
                  (const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi, VtxIndType* fv) const
{
    if ( (fi > CELL_DIM) ) std::cout << "ERROR: facet index fi is too small or too large!" << std::endl;
    assert( (fi <= CELL_DIM) );

    // note: vertex fi is opposite facet fi
    //       so, facet fi does NOT contain vertex fi
    for (SmallIndType kk=0; kk < fi; ++kk)
        fv[kk] = vtx[kk];
    for (SmallIndType kk=fi+1; kk < (CELL_DIM+1); ++kk)
        fv[kk-1] = vtx[kk];
}

/***************************************************************************************/
/* this returns the (facet) vertices in fv that are not equal to vi, where vi is also
   in the facet, i.e. it returns the vertices in the facet that are adjacent to vi.
   Note: if vi is not in the facet, then adj_vtx contains NULL_Vtx's (NULL values). */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Adj_Vertices_In_Facet
                  (const VtxIndType* fv, const VtxIndType& vi, VtxIndType* adj_vtx) const
{
    SmallIndType ii = NULL_Small; // init to invalid value
    for (SmallIndType kk=0; kk < CELL_DIM; ++kk)
        if (fv[kk]==vi)
        {
            ii = kk; // save it
            break;
        }

    // if we found a match
    if (ii != NULL_Small)
    {
        // copy array over (except for the matching entry)
        for (SmallIndType kk=0; kk<ii; ++kk)
            adj_vtx[kk] = fv[kk];
        for (SmallIndType kk=ii+1; kk<CELL_DIM; ++kk)
            adj_vtx[kk-1] = fv[kk];
    }
    else // return NULL value
    {
		for (int jj=0; jj < (int) (CELL_DIM-1); ++jj) // need (int) to avoid compiler warning!
			adj_vtx[jj] = NULL_Vtx;
    }
    // Note: if CELL_DIM==1, then adj_vtx has zero length, and nothing here will execute.
}

/***************************************************************************************/
/* Given the vertices of a cell and local facet index (indexing from 0),
   find the global vertex in that facet with the largest global index. */
template <SmallIndType CELL_DIM>
VtxIndType BM<CELL_DIM>::Get_Vertex_With_Largest_Index_In_Facet(const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi) const
{
    if ( (fi > CELL_DIM) ) std::cout << "ERROR: facet index fi is too small or too large!" << std::endl;
    assert( (fi <= CELL_DIM) );
    VtxIndType MAX = 0; // init to smallest possible value

    // loop through the vertices/facets
    // note: facet j (0 <= j < N) is opposite vertex j
    //       so, vertex j is NOT contained in facet j
    for (SmallIndType kk=0; kk < fi; ++kk)
        if (vtx[kk] > MAX) MAX = vtx[kk];
    for (SmallIndType kk=fi+1; kk < (CELL_DIM+1); ++kk)
        if (vtx[kk] > MAX) MAX = vtx[kk];

    return (const VtxIndType) MAX;
}

/***************************************************************************************/
/* return true if arrays are equal, otherwise false.
   Note: if the arrays have zero length, this returns true. */
template <SmallIndType CELL_DIM>
inline bool BM<CELL_DIM>::Adj_Vertices_In_Facet_Equal(const VtxIndType* a, const VtxIndType* b) const
{
    for (int ii=0; ii < (int) (CELL_DIM-1); ++ii)
        if (a[ii]!=b[ii]) return false;

    return true;
}

/* TODO */

// need to do a uniform hierarchical mesh refinement implementation

// need a multi-dim mesh manager, that can store sub-domains (cells), etc... and that can find sub-domain embeddings...

// have a class that manages multi-levels of BaseMesh...

// should have a routine to check the validity of the data structures...
// or some other sanity checks!!!!


// Here are the steps:
// 1. Finish developing the AHF mesh class.  *Everything*
// depends on that.  Also, if you give up on the FELICITY-
// Python conversion, at least you still have a product.
// Note: this includes adaptive mesh refinement, and having
// sub-domains.
// 2. Interface the mesh class to Python.  Again, this is a
// nice product to have by itself.
// 3. Then develop the DSL to only define domains, sub-domains,
// and some trivial function spaces, e.g. the Constant space


#undef BM

/***/
