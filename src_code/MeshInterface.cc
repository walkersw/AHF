/*
============================================================================================
   Interface (abstract) class for AHF Mesh class (of any dimensions).

Why did I do it like this?  Just have Mesh inherit from BaseMesh, and it will have an additional member variable "Points".

You did it so Mesh_Open could be used by BaseMesh and BasePtCoord.

But you can re-define Open in Mesh so that it affects the Mesh_Open sub-field of Points.

   Copyright (c) 12-17-2016,  Shawn W. Walker
============================================================================================
*/

#define _MESHINTERFACE_CC

#ifndef _PRELIM_H
#include "Prelim.h"  // simple typedefs, etc...
#endif

NOT USED!!!!

#ifdef _WIN32
/* 'class1' : inherits 'class2::member' via dominance */
#pragma warning( disable: 4250 ) // remove this annoying warning
#endif

/* C++ class definition */
#define  MMI  MeshInterface
class MMI
{
public:
    MMI();
    virtual ~MMI()=0;

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

    /* virtual methods for everything we want to access */
    virtual void Clear()=0;

    /* virtual methods for BaseMesh */

    // get the topological dimension
    virtual SmallIndType Top_Dim() const = 0;
    // allocate room for specified number of elements
    virtual void Reserve_Cells(const CellIndType&)=0;
    // append one cell to the end of the array
    virtual void Append_Cell(const VtxIndType*)=0;
    virtual void Append_Cell(const VtxIndType&, const VtxIndType&)=0;
    virtual void Append_Cell(const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    virtual void Append_Cell(const VtxIndType&, const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    // ... and update the intermediate data structure v2hfs (build incrementally)
    virtual void Append_Cell_And_Update(const VtxIndType*)=0;
    virtual void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&)=0;
    virtual void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    virtual void Append_Cell_And_Update(const VtxIndType&, const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    // get number of elements stored
    virtual CellIndType Num_Cells() const = 0;

    // retrieve one cell (writeable)
    virtual VtxIndType* Get_Cell_vtx(const CellIndType&)=0;
    virtual HalfFacetType* Get_Cell_halffacet(const CellIndType&)=0;
    virtual void Get_Cell(const CellIndType&, VtxIndType*, HalfFacetType*)=0;
    // retrieve one cell in its struct form (writeable)
    virtual CellType& Get_Cell(const CellIndType&)=0;

    // retrieve one cell (read-only)
    virtual const VtxIndType* Get_Cell_vtx(const CellIndType&) const = 0;
    virtual const HalfFacetType* Get_Cell_halffacet(const CellIndType&) const = 0;
    virtual void Get_Cell(const CellIndType&, const VtxIndType*&, const HalfFacetType*&) const = 0;
    // retrieve one cell in its struct form (read-only)
    virtual const CellType& Get_Cell(const CellIndType&) const = 0;

    // finalize the data structures for determining mesh connectivity
    //    (i.e. neighbors, vtx2half-facet mapping, etc.) and *close* the mesh.
    virtual void Finalize_Mesh_Connectivity()=0;

    /* all public routines below this need the mesh to be finalized to output
           correct information, i.e. the mesh should be "closed" and all internal
           data structures updated. This is done by building the sibling half-facet
           structure, and filling out the Vtx2HalfFacets mapping. All of this is
           automatically done by the "Finalize_Mesh_Connectivity" method. */

    // return const reference to v2hfs (an intermediate, internal data structure)
    virtual const Vtx2HalfFacet_Mapping& Get_v2hfs() const = 0;
    // get read access to Vtx2HalfFacets
    virtual const Vtx2HalfFacet_Mapping& Get_Vtx2HalfFacets() const = 0;

    // get unique set of vertices
    virtual void Get_Unique_Vertices(std::vector<VtxIndType>& uv) const = 0;
    virtual void Display_Unique_Vertices() const = 0;
    // get number of vertices referenced in Cell (don't forget to run Build_Vtx2HalfFacets)
    virtual VtxIndType Num_Vtx() const = 0;
    // get maximum vertex index referenced in Cell (don't forget to run Build_Vtx2HalfFacets)
    virtual VtxIndType Max_Vtx_Index() const = 0;

    // returns a unique set of all edges of the mesh
    virtual void Get_Edges(std::vector<MeshEdgeType>&) const = 0;
    // test if a pair of vertices is connected by an edge
    virtual bool Is_Connected(const VtxIndType&, const VtxIndType&) const = 0;
    virtual bool Is_Connected(const VtxIndType vv[2]) const = 0;
    // returns all cell indices attached to a given edge
    virtual void Get_Cells_Attached_To_Edge(const MeshEdgeType&, std::vector<CellIndType>&) const = 0;
    // returns all half-facets that are referenced by only one cell;
    //         i.e. the half-facets that are on the boundary of the mesh
    virtual void Get_FreeBoundary(std::vector<HalfFacetType>&) const = 0;

    // print out cell connectivity and sibling half-facet data
    virtual void Display_Cell(const CellIndType& ci=NULL_Cell) const = 0;
    // print out half-facets attached to vertex in intermediate data structure
    virtual void Display_v2hfs(const VtxIndType& vi=NULL_Vtx) const = 0;
    // print out half-facets attached to vertex for final data structure
    virtual void Display_Vtx2HalfFacets(const VtxIndType& vi=NULL_Vtx) const = 0;
    // returns all cell indices attached to a given vertex
    virtual void Get_Cells_Attached_To_Vertex(const VtxIndType&, std::vector<CellIndType>&) const = 0;
    // returns all cell indices attached to a given vertex and facet-connected to the given cell
    virtual void Get_Cells_Attached_To_Vertex(const VtxIndType&, const CellIndType&, std::vector<CellIndType>&) const = 0;
    // display cell indices attached to given vertex
    virtual void Display_Cells_Attached_To_Vertex(const VtxIndType&) const = 0;
    // determine if two cells (that share a vertex) are connected by a sequence of facet neighbors
    virtual bool Two_Cells_Are_Facet_Connected(const VtxIndType&, const CellIndType&, const CellIndType&) const = 0;
    // display if two cells (that both contain the given vertex) are facet-sequence-connected
    virtual void Display_Two_Cells_Are_Facet_Connected(const VtxIndType&, const CellIndType&, const CellIndType&) const = 0;
    // get all half-facets attached to a given half-facet
    virtual void Get_HalfFacets_Attached_To_HalfFacet(const HalfFacetType&, std::vector<HalfFacetType>&) const = 0;
    // display attached half-facets
    virtual void Display_HalfFacets_Attached_To_HalfFacet(const HalfFacetType&) const = 0;
    // get a unique set of non-manifold facets
    virtual void Get_Nonmanifold_HalfFacets(std::vector<HalfFacetType>&) const = 0;
    // display all non-manifold facets
    virtual void Display_Nonmanifold_HalfFacets() const = 0;
    // get the set of non-manifold vertices
    virtual void Get_Nonmanifold_Vertices(std::vector<VtxIndType>&) const = 0;
    // display all non-manifold vertices
    virtual void Display_Nonmanifold_Vertices() const = 0;

    /* these methods here are for debuggin purposes only.
       the casual user should never use them! */
    virtual void Finalize_v2hfs_DEBUG(bool bfs=true)=0;
    virtual void Build_Sibling_HalfFacets_DEBUG()=0;
    virtual void Build_Vtx2HalfFacets_DEBUG()=0;

    /* virtual methods for BasePtCoord */

    // get the topological dimension
    virtual SmallIndType Geo_Dim() const = 0;
    // get number of vertex coordinates stored
    virtual VtxIndType Num_Points() const = 0;
    // allocate room for specified number of vertex point coordinates (plus some room to grow)
    virtual void Reserve_Points(const VtxIndType&)=0;
    // initialize point coordinates to the *origin* for a specific number of vertices
    virtual void Init_Points(const VtxIndType&)=0;
    // set point coordinates of specific vertex
    virtual void Set_Coord(const VtxIndType& vi, const PointType* vc)=0;
    virtual void Set_Coord(const VtxIndType&, const PointType&)=0;
    virtual void Set_Coord(const VtxIndType&, const PointType&, const PointType&)=0;
    virtual void Set_Coord(const VtxIndType&, const PointType&, const PointType&, const PointType&)=0;

    // // retrieve one point (writeable)
    virtual PointType* Get_Point_coord(const VtxIndType&)=0;
    virtual CoordType& Get_Point(const VtxIndType&)=0;
    // retrieve one point (read-only)
    virtual const PointType* Get_Point_coord(const VtxIndType&) const = 0;
    virtual const CoordType& Get_Point(const VtxIndType&) const = 0;

    // print out vertex coordinates
    virtual void Display_Vtx_Coord(const VtxIndType& vi=NULL_Vtx) const = 0;

protected:
    /* flag to indicate if mesh cells and points may be added.
       true  = cells can be added, modified; point coordinate may be changed, etc.
       false = the mesh cells and points cannot be changed! */
    bool Mesh_Open;
};

/***************************************************************************************/
/* constructor */
MMI::MMI()
{
    //std::cout << "MeshInterface constructor..." << std::endl;
    Mesh_Open = true; // the mesh starts out as open for modification
}

/***************************************************************************************/
/* DE-structor */
MMI::~MMI()
{
    //std::cout << "MeshInterface destructor..." << std::endl;
}

// SWW: the only things that go here are interface functions to
//      BaseMesh and BasePtCoord

#undef MMI

/***/
