/*
============================================================================================
   Base class for array based half-facet (AHF) data structure to store and process meshes.

   Note: this base class is used in deriving the 1-D, 2-D, and 3-D mesh classes,
   as well as higher dimensions.
   Note: no vertex coordinates are stored in this class or the derived classes.

   Copyright (c) 11-05-2015,  Shawn W. Walker
============================================================================================
*/

#include "Prelim.h" // basic typedefs and includes
#include "Vtx2HalfFacet_Mapping.cc" // class for handling the Vertex-to-Half-Facet Mapping

/***************************************************************************************/
/* cell (simplex) connectivity data and sibling half-facets.
   CELL_DIM is the topological dimension of the cell (which is assumed to be a simplex).
   CELL_DIM = 1: cell is an edge (line segment);
                 half-facets are end-point vertices of the edge.
   CELL_DIM = 2: cell is a triangle;
                 half-facets are edges of the triangle.
   CELL_DIM = 3: cell is a tetrahedron;
                 half-facets are faces of the tetrahedron.
*/
template<SmallIndType CELL_DIM>
struct CellSimplexType
{
    VtxIndType           vtx[CELL_DIM+1]; // global vertex indices of the cell
    HalfFacetType  halffacet[CELL_DIM+1]; // half-facets corresponding to local facets of cell
                                          // note: see examples in concrete sub-class
    // check equality
    inline bool Equal(const CellSimplexType<CELL_DIM>& IN) const
    {
        for (SmallIndType ii = 0; ii < (CELL_DIM+1); ii++)
        {
            if (IN.vtx[ii]!=vtx[ii]) return false;
            if (!IN.halffacet[ii].Equal(halffacet[ii])) return false;
        }
        return true;
    }
    // set all values
    inline void Set(const SmallIndType& index, const VtxIndType& v, const HalfFacetType& hf)
    {
        vtx[index] = v;
        halffacet[index].Set(hf);
    }
};

/***************************************************************************************/
/* We need to define a one-to-one mapping between *local* facets and
   *local* vertices of a cell.

   For an edge (CELL_DIM==1), we define this to be:
      Vtx | Facet (Vertex)
     -----+-------------
       1  |   1
       2  |   2

   For a triangle (CELL_DIM==2), we define this to be:
      Vtx | Facet (Edge)
     -----+-------------
       1  |   2
       2  |   3
       3  |   1

   For a tetrahedron (CELL_DIM==3), we define this to be:
      Vtx | Facet (Face)
     -----+-------------
       1  |   2
       2  |   3
       3  |   4
       4  |   1

   (The pattern is obvious for higher cell dimensions...)

   Note: each vertex is contained (attached) to its corresponding facet.
*/
/* struct for "marking" which vertices/half-facets have been visited. */
template<SmallIndType CELL_DIM>
struct CellSimplexMarkedType
{
    bool      facet[CELL_DIM+1]; // local facets (true/false = visited/or not)
};

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
    void Clear() // clear all data
    {
        Cell.clear();
        Vtx2HalfFacets.Clear();
        v2hfs.Clear();
    };
    // allocate room for specified number of elements
    void Reserve_Cell(const CellIndType&);
    // append one cell to the end of the array
    void Append_Cell_Only(const VtxIndType v[(CELL_DIM+1)]);
    void Append_Cell_Only(const VtxIndType&, const VtxIndType&, const VtxIndType&);
    // ... and update the intermediate data structure v2hfs (build incrementally)
    void Append_Cell(const VtxIndType v[(CELL_DIM+1)]);
    void Append_Cell(const VtxIndType&, const VtxIndType&, const VtxIndType&);
    // get number of elements stored
    CellIndType Num_Cell() const { return (CellIndType) Cell.size(); };
    // retrieve one cell
    inline CellSimplexType<CELL_DIM>& Get_Cell(const CellIndType&); // for writeable
    inline const CellSimplexType<CELL_DIM>& Get_Cell(const CellIndType&) const; // for read-only
    // finalize intermediate data structure
    // (input = true means build from scratch using cell connectivity data)
    void Finalize_v2hfs(bool bfs=false);
    // return const reference to v2hfs (internal data)
    const Vtx2HalfFacet_Mapping& Get_v2hfs() const { const Vtx2HalfFacet_Mapping& c_v2hfs = v2hfs; return c_v2hfs; };
    // finishing filling out the sibling half-facet data structure
    void Build_Sibling_HalfFacets();
    // build the final Vtx2HalfFacets data struct
    void Build_Vtx2HalfFacets();

    /* all public routines below this need sibling half-facets and Vtx2HalfFacets
       to be completely finalized before they can be used! */
    // get unique set of vertices
    void Get_Unique_Vertices(std::vector<VtxIndType>& uv) const { Vtx2HalfFacets.Get_Unique_Vertices(uv); };
    void Display_Unique_Vertices() const { Vtx2HalfFacets.Display_Unique_Vertices(); };
    // get number of vertices referenced in Cell
    VtxIndType Num_Vtx() const
    {
        std::vector<VtxIndType> uv;
        Get_Unique_Vertices(uv);
        return (VtxIndType) uv.size();
    };
    // print out cell connectivity and sibling half-facet data
    void Display_Cell(const CellIndType& ci=NULL_IND) const;
    // print out half-facets attached to vertex in intermediate data structure
    void Display_v2hfs(const VtxIndType& vi=NULL_IND) const;
    // print out half-facets attached to vertex for final data structure
    void Display_Vtx2HalfFacets(const VtxIndType& vi=NULL_IND) const;
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
    // get a unique set of non-manifold facets
    void Get_Nonmanifold_HalfFacets(std::vector<HalfFacetType>&) const;
    // display all non-manifold facets
    void Display_Nonmanifold_HalfFacets() const;
    // get the set of non-manifold vertices
    void Get_Nonmanifold_Vertices(std::vector<VtxIndType>&) const;
    // display all non-manifold vertices
    void Display_Nonmanifold_Vertices() const;

    // connectivity and sibling half-facet data
    typedef CellSimplexType<CELL_DIM> CellSimplex_DIM;
    std::vector<CellSimplex_DIM>  Cell;
    // referenced vertices in Cell and (possibly multiple) attached half-facet(s)
    Vtx2HalfFacet_Mapping         Vtx2HalfFacets;

protected:
    double       Reserve_Buffer; // amount of extra memory to allocate when re-allocating
                                 // Cell and Vtx2HalfFacets (number between 0.0 and 1.0).
    VtxIndType   Estimate_Size_Vtx2HalfFacets; // estimate of the size to allocate in Vtx2HalfFacets.

    // intermediate data structure for building sibling half-facet information
    Vtx2HalfFacet_Mapping   v2hfs; // Note: for a given vertex, it references multiple half-facets.
    // Note: this data structure will NOT NECESSARILY store all referenced vertices
    //       in the triangulation.  This is because the vertex with smallest index
    //       will never be referenced (for example).  This is an internal structure that
    //       is only used to construct the sibling half-facet information (stored in Cell).

    // append half facets of given cell to v2hfs
    void Append_Half_Facets(const CellIndType&, const VtxIndType v[(CELL_DIM+1)]);
    // mapping used for generating the sibling half-facets
    void Vtx2Adjacent(const VtxIndType&, const CellIndType&, const SmallIndType&, VtxIndType va[CELL_DIM-1]) const;

    // amount extra to reserve when finding a variable number of cells attached to a vertex
    static const unsigned int cell_attach_chunk = 5*CELL_DIM;
    // recursion routines for public methods above
    void Get_Cells_Attached_To_Vertex_Recurse(const VtxIndType&, const CellIndType&, std::vector<CellIndType>&) const;
    // recursive call for public method above
    bool Two_Cells_Are_Facet_Connected_Recurse(const VtxIndType&, const CellIndType&, const CellIndType&,
                                               const CellIndType&, CellIndType&) const;

    /* baby routines */
    // given a global vertex, this determines which local vertex it is in the cell
    inline SmallIndType Get_Local_Vertex_Index_In_Cell(const VtxIndType&, const CellSimplexType<CELL_DIM>&) const;
    // map from local vertex in cell to local facets in cell that contain that vertex
    inline void Get_Local_Facets_Sharing_Local_Vertex(const SmallIndType& vi, SmallIndType facet[CELL_DIM]) const;
    // map local facet index to local vertices contained in that facet
    inline void Get_Local_Vertices_Of_Local_Facet(const SmallIndType&, SmallIndType v[CELL_DIM]) const;

    // get the global vertices in a cell's facet
    inline void Get_Global_Vertices_In_Facet
          (const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi, VtxIndType fv[CELL_DIM]) const;
    // get the vertices in a facet that are adjacent to the given vertex
    inline void Get_Adj_Vertices_In_Facet
          (const VtxIndType fv[CELL_DIM], const VtxIndType& vi, VtxIndType adj_vtx[(CELL_DIM-1)]) const;
    // simple max operation
    VtxIndType Get_Vertex_With_Largest_Index_In_Facet(const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi) const;
    // determine if adjacent vertices are equal
    inline bool Adj_Vertices_In_Facet_Equal(const VtxIndType a[CELL_DIM-1], const VtxIndType b[CELL_DIM-1]) const;

};

/***************************************************************************************/
/* constructor */
template <SmallIndType CELL_DIM>
BM<CELL_DIM>::BM ()
{
    Reserve_Buffer = 0.2; // allocate an extra 20% when re-allocating
    Estimate_Size_Vtx2HalfFacets = 0;
    // ensure memory is clear to start
    Clear();
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType CELL_DIM>
BM<CELL_DIM>::~BM()
{
    // clear the data
    Clear();
}

/***************************************************************************************/
/* Allocate memory to hold triangulation of given size (plus a little). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Reserve_Cell(const CellIndType& Num_Cells)
{
    // compute the actual size to allocate for the cells
    const CellIndType Actual_Cell_SIZE = (CellIndType) ((1.0 + Reserve_Buffer) * Num_Cells);
    Cell.reserve(Actual_Cell_SIZE);
    // guess on what to reserve for the intermediate data structure
    v2hfs.Reserve((VtxIndType) (CELL_DIM+1) * Num_Cells);
}

/***************************************************************************************/
/* Append a cell element to the end of the list.
   Note: input is the global vertex indices of the corners of the cell simplex. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_Only(const VtxIndType vtx[(CELL_DIM+1)])
{
    CellSimplexType<CELL_DIM>  CL;
    for (SmallIndType ii = 0; ii < (CELL_DIM+1); ii++)
    {
        CL.vtx[ii] = vtx[ii];
        // init to null value
        CL.halffacet[ii].Set();
    }

    /* append cell to the end */
    // first check that there is room; if not, then reserve more space
    const CellIndType Size_Est = (CellIndType) Cell.size();
    if (Size_Est >= (CellIndType) Cell.capacity())
        Reserve_Cell(Size_Est+2); // then out of room, so reserve more space
    Cell.push_back(CL);
}
/* Append a cell element to the end of the list (different format for input). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell_Only(const VtxIndType& v1, const VtxIndType& v2, const VtxIndType& v3)
{
    const VtxIndType vtx[CELL_DIM+1] = {v1, v2, v3};
    Append_Cell_Only(vtx);

}
/* Append a cell element to the end of the list,
   and build the intermediate v2hfs structure (incrementally). */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType vtx[(CELL_DIM+1)])
{
    Append_Cell_Only(vtx);

    // get current cell index
    const CellIndType ci = Cell.size(); // i.e. the current size
    Append_Half_Facets(ci, vtx);
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Append_Cell(const VtxIndType& v1, const VtxIndType& v2, const VtxIndType& v3)
{
    const VtxIndType vtx[CELL_DIM+1] = {v1, v2, v3};
    Append_Cell(vtx);
}

/***************************************************************************************/
/* store adjacent half-facets to vertices in intermediate data structure. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Finalize_v2hfs(bool Build_From_Scratch)
{
    if (Build_From_Scratch)
    {
        v2hfs.Clear(); // start fresh
        // record all vertex indices (with duplicates)
        const CellIndType NC = Num_Cell();
        v2hfs.Reserve((VtxIndType) (CELL_DIM+1) * NC);
        for (CellIndType ci = 1; ci <= NC; ci++)
        {
            const CellSimplexType<CELL_DIM>& CL = Cell[ci-1];
            Append_Half_Facets(ci, CL.vtx);
        }
    }
    // don't forget to sort!
    v2hfs.Sort();
}

/***************************************************************************************/
/* fill in the sibling half-facet data structure. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Build_Sibling_HalfFacets()
{
    // go thru all the elements
    const CellIndType NC = Num_Cell();
    for (CellIndType ci = 1; ci <= NC; ci++)
    {
        const CellIndType ci_C_style = ci - 1;
        // loop through all the facets
        for (SmallIndType ff = 1; ff <= (CELL_DIM+1); ff++)
        {
            if (Cell[ci_C_style].halffacet[ff-1].Is_Null()) // i.e. it is uninitialized
            {
                // find the vertex with largest ID in the face ff
                const VtxIndType MaxVtx = Get_Vertex_With_Largest_Index_In_Facet(Cell[ci_C_style].vtx, ff);
                // get vertices in the facet ff of the current element that is adjacent to MaxVtx
                VtxIndType  Adj_Vtx[CELL_DIM-1];
                Vtx2Adjacent(MaxVtx, ci, ff, Adj_Vtx); // see below...

                // find all half-facets that are attached to MaxVtx
                std::pair<std::vector<VtxHalfFacetType>::const_iterator,
                          std::vector<VtxHalfFacetType>::const_iterator>  RR; // make range variable
                const unsigned int Num_HF = v2hfs.Get_Half_Facets(MaxVtx, RR);

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
                        VtxIndType  Adj_Vtx_First[CELL_DIM-1];
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
                        VtxIndType  Adj_Vtx_Other[CELL_DIM-1];
                        Vtx2Adjacent(MaxVtx, (*Next).ci, (*Next).fi, Adj_Vtx_Other);

                        // if the half-facet is valid
                        if (Adj_Vertices_In_Facet_Equal(Adj_Vtx_Other, Adj_Vtx))
                        {
                            // in the Current cell and half-facet, write the Next half-facet
                            CellSimplexType<CELL_DIM>& Current_Cell = Get_Cell((*Current).ci);
                            Current_Cell.halffacet[(*Current).fi-1].ci = (*Next).ci;
                            Current_Cell.halffacet[(*Current).fi-1].fi = (*Next).fi;
                            // Update Current to Next
                            Current = Next;
                        }
                    }
                    // don't forget to close the cycle:
                    // if Current is different from Start
                    if (Current!=Start) // i.e. it cannot refer to itself!
                    {
                        // then in the Current cell and half-facet, write the Start half-facet
                        CellSimplexType<CELL_DIM>& Current_Cell = Get_Cell((*Current).ci);
                        Current_Cell.halffacet[(*Current).fi-1].ci = (*Start).ci;
                        Current_Cell.halffacet[(*Current).fi-1].fi = (*Start).fi;
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
    Vtx2HalfFacets.Clear(); // start fresh
    // make guess on how much space to allocate
    Vtx2HalfFacets.Reserve(Estimate_Size_Vtx2HalfFacets);

    // allocate a temp variable and initialize to all false,
    //    because we have not visited any half-facets yet.
	typedef CellSimplexMarkedType<CELL_DIM> CellMarked_DIM;
    CellMarked_DIM mm;
    for (SmallIndType kk = 0; kk < (CELL_DIM+1); ++kk)
        mm.facet[kk] = false;
    const CellIndType NC = Num_Cell();
    std::vector<CellMarked_DIM> Marked(NC,mm);

    // store one-to-one mapping from local vertices to local facets
    SmallIndType lv_to_lf[CELL_DIM+1];
    if (CELL_DIM==1)
    {
        lv_to_lf[0] = 1;
        lv_to_lf[1] = 2;
    }
    else
    {
        for (SmallIndType kk = 1; kk < (CELL_DIM+1); ++kk)
            lv_to_lf[kk-1] = kk+1;
        lv_to_lf[CELL_DIM] = 1;
    }

    // loop through each cell
    for (CellIndType ci = 1; ci <= NC; ci++)
    {
        // loop through each local vertex of the current cell
        for (SmallIndType local_vi = 1; local_vi <= (CELL_DIM+1); local_vi++)
        {
            // current half-facet is <ci,local_fi>, local_fi = lv_to_lf(local_vi)
            const SmallIndType local_fi = lv_to_lf[local_vi-1];
            if (!Marked[ci-1].facet[local_fi-1]) // need C-style indexing!
            {
                // store this half-facet
                VtxHalfFacetType vhf;
                const VtxIndType Global_Vtx = Cell[ci-1].vtx[local_vi-1]; // get global vertex
                vhf.vtx = Global_Vtx;
                vhf.ci  = ci;
                vhf.fi  = local_fi;
                Vtx2HalfFacets.Append(vhf);
                // denote that we have visited it!
                Marked[ci-1].facet[local_fi-1] = true;

                // get the cells attached to the current vertex, that are also
                //     facet-connected to the current cell
                std::vector<CellIndType> attached_cells;
                Get_Cells_Attached_To_Vertex(Global_Vtx, ci, attached_cells);
                // loop thru the attached cells and mark corresponding half-facets as also visited
                for (std::vector<CellIndType>::iterator it = attached_cells.begin(); it < attached_cells.end(); ++it)
                {
                    const CellIndType ci_hat = *it;
                    // the local vertex index within ci_hat
                    const SmallIndType local_vi_hat = Get_Local_Vertex_Index_In_Cell(Global_Vtx, Cell[ci_hat-1]);
                    // corresponding local facet index
                    const SmallIndType local_fi_hat = lv_to_lf[local_vi_hat-1];
                    // mark it!
                    Marked[ci_hat-1].facet[local_fi_hat-1] = true;
                }
            }
        }
    }

    Vtx2HalfFacets.Sort(); // now this data structure is usable!

    /* Give border half-facets (i.e. half-facets with no siblings) higher priority.
       This allows us to easily identify vertices that are on the boundary of the mesh. */

    // loop through each cell
    for (CellIndType ci = 1; ci <= NC; ci++)
    {
        // loop through each local facet of the current cell
        for (SmallIndType local_fi = 1; local_fi <= (CELL_DIM+1); local_fi++)
        {
            // if this half-facet has no sibling
            if (Cell[ci-1].halffacet[local_fi-1].Is_Null())
            {
                // get the local vertices of the local facet
                SmallIndType local_vtx[CELL_DIM];
                Get_Local_Vertices_Of_Local_Facet(local_fi, local_vtx);

                // for each vertex of the current half-facet
                for (SmallIndType local_vi = 1; local_vi <= CELL_DIM; local_vi++)
                {
                    // get the global vertex
                    const VtxIndType global_vi = Cell[ci-1].vtx[local_vtx[local_vi-1] - 1];

                    // get the half-facets attached to the vertex
                    std::pair <std::vector<VtxHalfFacetType>::iterator,
                               std::vector<VtxHalfFacetType>::iterator> RR;
                    const unsigned int Num_HF = Vtx2HalfFacets.Get_Half_Facets(global_vi, RR);

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
            }
        }
    }
    // Note: we don't have to sort again,
    //       because the half-facets are ordered by the attached vertex
}

/***************************************************************************************/
/* print cell connectivity and sibling haf-facets. "ci" is the index of a
   specific cell; if ci=NULL_IND, then print all cells. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Cell(const CellIndType& ci) const
{
    if (ci==NULL_IND)
    {
        // then print all the cells
        std::cout << "Display connectivity of all cells:" << std::endl;
        std::cout << "Cell #    |        Vertices        |   Sibling Half-Facets" << std::endl;
        for (CellIndType ii = 1; ii <= Num_Cell(); ii++)
        {
            const CellSimplexType<CELL_DIM>& CL = Get_Cell(ii);
            std::cout << ii << "  |  " << CL.vtx[0]; // print cell # and first vertex
            for (SmallIndType kk = 1; kk < (CELL_DIM+1); kk++)
                std::cout << ", " << CL.vtx[kk];
            std::cout << "   |   " << "<" << CL.halffacet[0].ci << "," << CL.halffacet[0].fi << ">";
            for (SmallIndType kk = 1; kk < (CELL_DIM+1); kk++)
                std::cout << ", " << "<" << CL.halffacet[kk].ci << "," << CL.halffacet[kk].fi << ">";
            std::cout << std::endl;
        }
    }
    else
    {
        // then print all the cells
        std::cout << "Display connectivity of cell #" << ci << ":" << std::endl;
        std::cout << "        Vertices        |   Sibling Half-Facets" << std::endl;

        const CellSimplexType<CELL_DIM>& CL = Get_Cell(ci);
        std::cout << "  " << CL.vtx[0]; // print first vertex
        for (SmallIndType kk = 1; kk < (CELL_DIM+1); kk++)
            std::cout << ", " << CL.vtx[kk];
        std::cout << "   |   " << "<" << CL.halffacet[0].ci << "," << CL.halffacet[0].fi << ">";
        for (SmallIndType kk = 1; kk < (CELL_DIM+1); kk++)
            std::cout << ", " << "<" << CL.halffacet[kk].ci << "," << CL.halffacet[kk].fi << ">";
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
/* Get reference to specific mesh cell data (given the cell index). */
template <SmallIndType CELL_DIM>
inline CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell(const CellIndType& ci)
{
    // ci must be between 1 and Num_Cell
    assert((ci>=1) && (ci<=Num_Cell()));

    const CellIndType Index = (ci - 1);
    return Cell[Index];
}
template <SmallIndType CELL_DIM>
inline const CellSimplexType<CELL_DIM>& BM<CELL_DIM>::Get_Cell(const CellIndType& ci) const // read-only
{
    // ci must be between 1 and Num_Cell
    assert((ci>=1) && (ci<=Num_Cell()));

    const CellIndType Index = (ci - 1);
    return Cell[Index];
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
    if (vi==NULL_IND) return; // the vertex is null, so do nothing.

    // get the attached half-facets
    std::pair <std::vector<VtxHalfFacetType>::const_iterator,
               std::vector<VtxHalfFacetType>::const_iterator> RR;
    const unsigned int Num_HF = Vtx2HalfFacets.Get_Half_Facets(vi, RR);

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
        cell_array.resize(old_size + temp_array.size(), NULL_IND);
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
    if ( (vi==NULL_IND) || (ci==NULL_IND) ) return; // the vertex or starting cell is null, so do nothing.

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
    if (ci==NULL_IND) return; // this can happen if the neighbor cell does not exist

    std::vector<CellIndType>::iterator it;
    it = std::find(cell_array.begin(), cell_array.end(), ci); // see if ci is already in the array
    // if ci is already in the array
    if (it!=cell_array.end())
        return; // then we have already visited this cell, so done
    else
    {
        // access the cell
        const CellSimplexType<CELL_DIM>& CL = Get_Cell(ci);
        // check again that the vertex is actually in the cell
        const SmallIndType local_vi = Get_Local_Vertex_Index_In_Cell(vi, CL);
        if (local_vi==NULL_IND) return; // cell does not actually contain the vertex (this should not happen)

        // add the cell to the list
        const CellIndType current_size = (CellIndType) cell_array.size();
        if (current_size == cell_array.capacity())
            cell_array.reserve(current_size + cell_attach_chunk);
        cell_array.push_back(ci);

        // get the local facets that share the vertex
        SmallIndType local_facet[CELL_DIM];
        Get_Local_Facets_Sharing_Local_Vertex(local_vi, local_facet);

        // loop through the facets
        for (SmallIndType fi = 0; fi < CELL_DIM; ++fi)
            if (local_facet[fi]!=NULL_IND)
            {
                // determine the neighbor cell on the "other side" of that facet and search it
                Get_Cells_Attached_To_Vertex_Recurse(vi, CL.halffacet[local_facet[fi]-1].ci, cell_array);
            }
            // else the neighbor does not exist, so do nothing
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
    if (vi==NULL_IND)
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
    if ( (vi==NULL_IND) || (ci_a==NULL_IND) || (ci_b==NULL_IND) ) return false; // something is null, so do nothing.

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
    if (current==NULL_IND) return false; // this can happen if the neighbor cell does not exist

    // if current matches the target
    if (current==target)
        return true; // then we found it!
    else
    {
        // access the cell
        const CellSimplexType<CELL_DIM>& CL = Get_Cell(current);
        // check again that the vertex is actually in the cell
        const SmallIndType local_vi = Get_Local_Vertex_Index_In_Cell(vi, CL);
        if (local_vi==NULL_IND) return false; // cell does not actually contain the vertex (this should not happen)

        // get the local facets that share the vertex
        SmallIndType local_facet[CELL_DIM];
        Get_Local_Facets_Sharing_Local_Vertex(local_vi, local_facet);

        // keep counting
        Depth_Count++;

        // loop through the facets
        for (SmallIndType fi = 0; fi < CELL_DIM; ++fi)
            if (local_facet[fi]!=NULL_IND)
            {
                // determine the neighbor cell on the "other side" of that facet and search it
                const bool CONNECTED = Two_Cells_Are_Facet_Connected_Recurse(vi, start, CL.halffacet[local_facet[fi]-1].ci, target, Depth_Count);
                if (CONNECTED) return true;
            }
            // else the neighbor does not exist, so do nothing

        return false; // else it was not found!
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
    if (vi==NULL_IND)
    {
        std::cout << "Vertex is invalid, so display nothing!" << std::endl;
        return; // the vertex is null, so do nothing.
    }
    if ( (ci_a==NULL_IND) || (ci_b==NULL_IND) )
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
        const CellSimplexType<CELL_DIM> CL_a = Get_Cell(ci_a);
        for (SmallIndType kk = 0; kk < (CELL_DIM+1); ++kk)
            if (CL_a.vtx[kk]==vi)
            {
                contain_a = true;
                break;
            }
        bool contain_b = false;
        const CellSimplexType<CELL_DIM> CL_b = Get_Cell(ci_b);
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
/* get all half-facets attached to a given half-facet.  This returns the list in a
   std::vector.  Note: the output also contains the given half-facet.
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
    unsigned int COUNT = 0;
    while (COUNT < 100000) // allow for up to 100,000 neighbors!
    {
        ++COUNT;
        // get the next half-facet
        const HalfFacetType& current_hf = attached_hf.back();
        const CellSimplexType<CELL_DIM> CL = Get_Cell(current_hf.ci);
        const HalfFacetType& next_hf = CL.halffacet[current_hf.fi-1];

        // if the neighbor does not exist, then stop!
        if (next_hf.Is_Null()) break;

        // if we get back to the starting half-facet, stop!
        if (next_hf.Equal(hf_in)) break;
        // else store it
        const unsigned int current_capacity = attached_hf.capacity();
        // reserve in chunks
        if (attached_hf.size()==current_capacity)
            attached_hf.reserve(current_capacity + 5);
        attached_hf.push_back(next_hf);
    }

    if (COUNT > 100000)
    {
        // then quit!
        std::cout << "ERROR in 'Get_HalfFacets_Attached_To_HalfFacet'..." << std::endl;
        std::cout << "    Number of neighbors is too great." << std::endl;
        std::cout << "    There should not be more than 100,000 cells attached to a single facet!" << std::endl;
    }
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_HalfFacets_Attached_To_HalfFacet(const HalfFacetType& hf_in) const
{
    std::vector<HalfFacetType> attached;
    Get_HalfFacets_Attached_To_HalfFacet(hf_in, attached);

    std::cout << "The half-facets attached to <" << hf_in.ci << "," << hf_in.fi << "> are:" << std::endl;
    for (unsigned int jj = 0; jj < attached.size(); ++jj)
        std::cout << "<" << attached[jj].ci << "," << attached[jj].fi << ">" << std::endl;
}

/***************************************************************************************/
/* get a unique set of non-manifold facets. Note: this requires the sibling
   half-facet data to be completed.  This returns a vector of half-facets, each
   defining a *distinct* non-manifold facet. */
bool non_manifold_sort_function(HalfFacetType A, HalfFacetType B) { return (A.ci < B.ci); };
bool non_manifold_equal_function(HalfFacetType& A, HalfFacetType& B) { return ( A.Equal(B) ); };
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Get_Nonmanifold_HalfFacets(std::vector<HalfFacetType>& non_manifold_hf) const
{
    non_manifold_hf.clear(); // start fresh
    non_manifold_hf.reserve(10); // init

    // go thru all the elements
    const CellIndType NC = Num_Cell();
    for (CellIndType ci = 1; ci <= NC; ++ci)
    {
        const CellSimplexType<CELL_DIM> CL = Get_Cell(ci);
        // loop through all the (half) facets
        for (SmallIndType fi = 1; fi <= (CELL_DIM+1); ++fi)
        {
            // get the facet neighbor
            const HalfFacetType& neighbor_hf = CL.halffacet[fi-1];
            const SmallIndType& n_ci = neighbor_hf.ci;
            const SmallIndType& n_fi = neighbor_hf.fi;

            // if the neighbor is not-NULL
            if (!neighbor_hf.Is_Null())
            {
                const CellSimplexType<CELL_DIM> N_CL = Get_Cell(n_ci);
                // if the neighbor half-facet looks back at the cell we started at
                if (N_CL.halffacet[n_fi-1].ci==ci)
                {
                    // and if the local facet does *not* match where we started
                    if (N_CL.halffacet[n_fi-1].fi!=fi)
                    {
                        // then: two distinct facets of a cell are joined together!!
                        // this should not happen!
                        std::cout << "ERROR in 'Get_Nonmanifold_HalfFacets':" << std::endl;
                        std::cout << "      Two facets of a cell are siblings; this should not happen!" << std::endl;
                        assert(N_CL.halffacet[n_fi-1].fi==fi);
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
                    const unsigned int current_capacity = non_manifold_hf.capacity();
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
    non_manifold_hf.resize( (unsigned int) std::distance(non_manifold_hf.begin(), it) );
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Nonmanifold_HalfFacets() const
{
    std::vector<HalfFacetType> non_manifold_hf;
    Get_Nonmanifold_HalfFacets(non_manifold_hf);

	const unsigned int NUM = non_manifold_hf.size();
	if (NUM==0)
		std::cout << "There are *no* non-manifold facets." << std::endl;
	else // there is at least 1
	{
		std::cout << "These are all the non-manifold facets in the mesh (output: <cell index, local facet index>):" << std::endl;
		for (unsigned int jj = 0; jj < NUM; ++jj)
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
            const unsigned int current_capacity = non_manifold_vtx.capacity();
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
    non_manifold_vtx.resize( (unsigned int) std::distance(non_manifold_vtx.begin(), stop) );
}
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Display_Nonmanifold_Vertices() const
{
    std::vector<VtxIndType> non_manifold_vtx;
    Get_Nonmanifold_Vertices(non_manifold_vtx);

	const unsigned int NUM = non_manifold_vtx.size();
	if (NUM==0)
		std::cout << "There are *no* non-manifold vertices." << std::endl;
	else // there is more than 1
	{
		std::cout << "These are all the non-manifold vertex indices:" << std::endl;
		for (unsigned int jj = 0; jj < NUM; ++jj)
			std::cout << non_manifold_vtx[jj] << std::endl;
	}
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
    for (SmallIndType fi = 1; fi <= (CELL_DIM+1); ++fi)
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
   If the half-facet does not contain the given vertex, then the output is 0.
   Note: <ci,fi> is a half-facet, where ci is a cell index, and fi is a local
   facet index. */
template <SmallIndType CELL_DIM>
void BM<CELL_DIM>::Vtx2Adjacent(const VtxIndType& v_in, const CellIndType& ci, const SmallIndType& fi,
                                VtxIndType v_adj[CELL_DIM-1]) const
{
    if ( (fi < 1) || (fi > (CELL_DIM+1)) ) std::cout << "ERROR: local facet index is invalid!" << std::endl;
    assert( (fi >= 1) && (fi <= (CELL_DIM+1)) );

    if (CELL_DIM>=2)
    {
        // get the vertices in the half-facet "fi"
        VtxIndType vtx_in_hf[CELL_DIM];
        const CellSimplexType<CELL_DIM> CL = Get_Cell(ci);
        Get_Global_Vertices_In_Facet(CL.vtx, fi, vtx_in_hf);

        // now return the vertices in the half-facet, EXCEPT v_in (i.e. the adjacent vertices)
        Get_Adj_Vertices_In_Facet(vtx_in_hf, v_in, v_adj);
    }
    // else do nothing!
}

/***************************************************************************************/
/* Find the *local* index of the given *global* vertex within the given cell. */
template <SmallIndType CELL_DIM>
inline SmallIndType BM<CELL_DIM>::Get_Local_Vertex_Index_In_Cell(const VtxIndType& vi, const CellSimplexType<CELL_DIM>& CL) const
{
    for (SmallIndType ii=1; ii<=(CELL_DIM+1); ++ii)
        if (CL.vtx[ii-1]==vi) return ii;

    return NULL_IND;
}

/***************************************************************************************/
/* return the local facet indices that share a given local vertex.
   Note: indexing starts at 1. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Local_Facets_Sharing_Local_Vertex(const SmallIndType& vi, SmallIndType facet[CELL_DIM]) const
{
    if ((vi < 1) || (vi > (CELL_DIM+1))) std::cout << "ERROR: vertex index vi is too small or too large!" << std::endl;
    assert((vi>=1)&&(vi<=(CELL_DIM+1)));

    // note: vertex vi is opposite facet vi
    //       so, facet vi does NOT contain vertex vi

    // get the local facets that contain local vertex vi
    for (SmallIndType kk=1; kk<vi; ++kk)
        facet[kk-1] = kk;
    for (SmallIndType kk=vi+1; kk<=(CELL_DIM+1); ++kk)
        facet[kk-2] = kk;
    // note: be careful with C-style indexing...
}

/***************************************************************************************/
/* Return the local vertices that are attached to a given local facet.
   Note: indexing starts at 1. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Local_Vertices_Of_Local_Facet(const SmallIndType& fi, SmallIndType vert[CELL_DIM]) const
{
    // hahah, we can reuse this!
    Get_Local_Facets_Sharing_Local_Vertex(fi, vert);
}

/***************************************************************************************/
/* given a cell's vertices and a local facet, return the global vertices contained
   in that facet. Note: fi indexes from 1. */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Global_Vertices_In_Facet
                  (const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi, VtxIndType fv[CELL_DIM]) const
{
    if ((fi < 1) || (fi > (CELL_DIM+1))) std::cout << "ERROR: facet index fi is too small or too large!" << std::endl;
    assert((fi>=1)&&(fi<=(CELL_DIM+1)));

    // note: vertex fi is opposite facet fi
    //       so, facet fi does NOT contain vertex fi
    for (SmallIndType kk=0; kk<(fi-1); ++kk)
        fv[kk] = vtx[kk];
    for (SmallIndType kk=fi; kk<(CELL_DIM+1); ++kk)
        fv[kk-1] = vtx[kk];
    // note: be careful with C-style indexing...
}

/***************************************************************************************/
/* this returns the (facet) vertices in fv that are not equal to vi, where vi is also
   in the facet, i.e. it returns the vertices in the facet that are adjacent to vi.
   Note: if vi is not in the facet, then adj_vtx contains NULL_IND's (NULL values). */
template <SmallIndType CELL_DIM>
inline void BM<CELL_DIM>::Get_Adj_Vertices_In_Facet
                  (const VtxIndType fv[CELL_DIM], const VtxIndType& vi, VtxIndType adj_vtx[(CELL_DIM-1)]) const
{
    SmallIndType ii = CELL_DIM; // init to invalid value
    for (SmallIndType kk=0; kk<CELL_DIM; ++kk)
        if (fv[kk]==vi)
        {
            ii = kk; // save it
            break;
        }

    // if we found a match
    if (ii < CELL_DIM)
    {
        // copy array over (except for the matching entry)
        for (SmallIndType kk=0; kk<ii; ++kk)
            adj_vtx[kk] = fv[kk];
        for (SmallIndType kk=ii+1; kk<CELL_DIM; ++kk)
            adj_vtx[kk-1] = fv[kk];
        // note: be careful with C-style indexing...
    }
    else // return NULL value
    {
        for (SmallIndType jj=0; jj<(CELL_DIM-1); ++jj)
            adj_vtx[jj] = NULL_IND;
    }
}

/***************************************************************************************/
/* Given the vertices of a cell and local facet index (indexing from 1),
   find the vertex in that facet with the largest index. */
template <SmallIndType CELL_DIM>
VtxIndType BM<CELL_DIM>::Get_Vertex_With_Largest_Index_In_Facet(const VtxIndType vtx[(CELL_DIM+1)], const SmallIndType& fi) const
{
    VtxIndType MAX = 0; // init to smallest possible value

    // loop through the vertices/facets
    // note: facet j (1 <= j <= N) is opposite vertex j
    //       so, vertex j is NOT contained in facet j
    for (SmallIndType kk=0; kk<(CELL_DIM+1); ++kk)
    {
        // ignore the vertex that is NOT in the facet
        if (kk!=(fi-1))
            if (vtx[kk] > MAX) MAX = vtx[kk];
    }
    return (const VtxIndType) MAX;
}

/***************************************************************************************/
/* return true if arrays are equal, otherwise false.
   Note: if the arrays have zero length, this returns true. */
template <SmallIndType CELL_DIM>
inline bool BM<CELL_DIM>::Adj_Vertices_In_Facet_Equal(const VtxIndType a[CELL_DIM-1], const VtxIndType b[CELL_DIM-1]) const
{
    for (SmallIndType ii=0; ii<(CELL_DIM-1); ++ii)
        if (a[ii]!=b[ii]) return false;

    return true;
}

/* TODO */

// neighbor query; this is just the sibling half-facets....

// store sub-domain cells...

// eventually, need a multi-dim mesh manager that can find sub-domain embeddings...

// put in functionality that TriRep has...

// should have a routine to check the validity of the data structures...
// or some other sanity checks!!!!


#undef BM

/***/
