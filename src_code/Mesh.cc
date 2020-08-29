/*
============================================================================================
   Class for array based half-facet (AHF) data structure.  This is a sub-class of BaseMesh
   which stores a *pointer* to a BasePtCoord object.  This way, many different meshes,
   of varying topological dimension, can share the same vertex (point) coordinates.  This
   also means that we can implement various methods in a more general way.

   Note: Everything is indexed starting at 0!

   Also note: using a vector of structs is 2 x faster than using a vector of integers.

   Copyright (c) 05-23-2020,  Shawn W. Walker
============================================================================================
*/

#define _MESH_CC

#ifndef _SIMPLEXMATH_H
#include "SimplexMath.h"  // math operations related to simplices
#endif
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
class MMC: public BaseMesh<CELL_DIM>
{
public:
    MMC();
    MMC(BasePtCoord<GEO_DIM>*);
    ~MMC();

    // open Mesh for modification
    inline void Open()
    {
        BaseMesh<CELL_DIM>::Open();
        Vtx->Open();
    };
    // close Mesh; modification is no longer allowed
    inline void Close()
    {
        BaseMesh<CELL_DIM>::Close();
        Vtx->Close();
    };

    // clear all data
    void Clear()
    {
        BaseMesh<CELL_DIM>::Clear();
        Vtx->Clear();
    };

    // set pointer to vertex coordinates
    void Set_Vertex_Data(BasePtCoord<GEO_DIM>* Input_Vtx) { Vtx = Input_Vtx; }
    // get pointer to vertex coordinates
    inline BasePtCoord<GEO_DIM>* Get_Vertex_Data_Ptr() { return (BasePtCoord<GEO_DIM>*) Vtx; };

    // coordinate conversion
    void Reference_To_Cartesian(const CellIndType&, const CellIndType*, const PointType*, PointType*);
    void Cartesian_To_Reference(const CellIndType&, const CellIndType*, const PointType*, PointType*);
    void Barycentric_To_Reference(const CellIndType&, const PointType*, PointType*);
    void Reference_To_Barycentric(const CellIndType&, const PointType*, PointType*);
    void Barycentric_To_Cartesian(const CellIndType&, const CellIndType*, const PointType*, PointType*);
    void Cartesian_To_Barycentric(const CellIndType&, const CellIndType*, const PointType*, PointType*);

    // cell quantities
    void Diameter(const CellIndType&, const CellIndType*, RealType*);
    void Bounding_Box(const CellIndType&, const CellIndType*, PointType*, PointType*);
    void Bounding_Box(PointType*, PointType*);
    void Volume(const CellIndType&, const CellIndType*, RealType*);
    void Angles(const CellIndType&, const CellIndType*, RealType*);

    // simplex centers
    void Barycenter(const CellIndType&, const CellIndType*, PointType*);
    void Circumcenter(const CellIndType&, const CellIndType*, PointType*, RealType*);
    void Incenter(const CellIndType&, const CellIndType*, PointType*, RealType*);
    
    // others
    void Shape_Regularity(const CellIndType&, const CellIndType*, RealType*);
    
private:

    /* keep a pointer to the set of global vertex coordinates */
    BasePtCoord<GEO_DIM>*  Vtx;

    // basic internal method
    void init();
};

/***************************************************************************************/
/* constructor */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
MMC<CELL_DIM, GEO_DIM>::MMC() : BaseMesh<CELL_DIM>()
{
    init(); // basic initialization
}

/***************************************************************************************/
/* constructor */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
MMC<CELL_DIM, GEO_DIM>::MMC(BasePtCoord<GEO_DIM>*  Input_Vtx) : BaseMesh<CELL_DIM>()
{
    init(); // basic initialization
    Set_Vertex_Data(Input_Vtx);
    Open(); // init to open
}

/***************************************************************************************/
/* set pointer for the vertex coordinate data */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::init()
{
    // initialize pointer to vertex coordinates
    Vtx = (BasePtCoord<GEO_DIM>*) NULL;

    // what else to do or check?
    if (GEO_DIM < CELL_DIM)
    {
        std::cout << "Desired topological dimension of a cell is " << CELL_DIM << "." << std::endl;
        std::cout << "Desired geometric dimension is " << GEO_DIM << "." << std::endl;
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
    Clear();
    //std::cout << "Mesh destructor..." << std::endl;
}

// /***************************************************************************************/
// /* set pointer for the vertex coordinate data */
// template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
// void MMC<CELL_DIM, GEO_DIM>::Set_Vertex_Data(BasePtCoord<GEO_DIM>*  Input_Vtx)
// {
    // // set pointer to vertex coordinates
    // Vtx = Input_Vtx;
// }

/***************************************************************************************/
/* convert reference element coordinates to cartesian coordinates, where the reference
   element is the "standard" reference simplex.
   Inputs: number of cells, and the cell indices (an array),
           reference coordinates (an array of consecutive groups of length CELL_DIM
   Outputs: cartesian coordinates (an array of consecutive groups of length GEO_DIM
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Reference_To_Cartesian(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             const PointType* PR,
                             PointType* PC)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        
        // get the relevant ref-coord for the current cell
        const PointType* PR_ii = PR + ii*CELL_DIM;
        // get relevant output for current cell
        PointType* PC_local = PC + ii*GEO_DIM;
        
        Simplex_Reference_To_Cartesian<CELL_DIM, GEO_DIM>(c_Vtx, VI, PR_ii, PC_local);
    }
}

/***************************************************************************************/
/* convert cartesian coordinates to reference element coordinates, where the reference
   element is the "standard" reference simplex.
   Inputs: number of cells, and the cell indices (an array),
           cartesian coordinates (an array of consecutive groups of length GEO_DIM
   Outputs: reference coordinates (an array of consecutive groups of length CELL_DIM
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Cartesian_To_Reference(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             const PointType* PC,
                             PointType* PR)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);

        // get relevant cart-coord for current cell
        const PointType* PC_local = PC + ii*GEO_DIM;
        // get the relevant output for the current cell
        PointType* PR_ii = PR + ii*CELL_DIM;
        
        Simplex_Cartesian_To_Reference<CELL_DIM, GEO_DIM>(c_Vtx, VI, PC_local, PR_ii);
    }
}

/***************************************************************************************/
/* convert barycentric coordinates to reference element coordinates.
   Inputs: number of cells, 
           barycentric coordinates (an array of consecutive groups of length (CELL_DIM+1)
   Outputs: reference coordinates (an array of consecutive groups of length CELL_DIM
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Barycentric_To_Reference(
                             const CellIndType& Num_Cell, const PointType* PB,
                             PointType* PR)
{
    // just copy over the 1st thru CELL_DIM barycentric coordinates
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get relevant bary-coord for current cell
        const PointType* PB_local = PB + ii*(CELL_DIM+1);
        // get the relevant output for the current cell
        PointType* PR_ii = PR + ii*CELL_DIM;
        
        Simplex_Barycentric_To_Reference<CELL_DIM>(PB_local, PR_ii);
    }
}

/***************************************************************************************/
/* convert reference element coordinates to barycentric coordinates.
   Inputs: number of cells, 
           reference coordinates (an array of consecutive groups of length CELL_DIM
   Outputs: barycentric coordinates (an array of consecutive groups of length (CELL_DIM+1)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Reference_To_Barycentric(
                             const CellIndType& Num_Cell, const PointType* PR,
                             PointType* PB)
{
    // just copy over the 1st thru (CELL_DIM) barycentric coordinates
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get relevant ref-coord for current cell
        const PointType* PR_local = PR + ii*CELL_DIM;
        // get the relevant output for the current cell
        PointType* PB_ii = PB + ii*(CELL_DIM+1);
        
        Simplex_Reference_To_Barycentric<CELL_DIM>(PR_local, PB_ii);
    }
}

/***************************************************************************************/
/* convert barycentric coordinates to cartesian coordinates.
   Inputs: number of cells, and the cell indices (an array),
           barycentric coordinates (an array of consecutive groups of length (CELL_DIM+1)
   Outputs: cartesian coordinates (an array of consecutive groups of length GEO_DIM
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Barycentric_To_Cartesian(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             const PointType* PB,
                             PointType* PC)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        
        // get the relevant bary-coord for the current cell
        const PointType* PB_local = PB + ii*(CELL_DIM+1);
        // get relevant output for current cell
        PointType* PC_ii = PC + ii*GEO_DIM;
        
        Simplex_Barycentric_To_Cartesian<CELL_DIM, GEO_DIM>(c_Vtx, VI, PB_local, PC_ii);
    }
}

/***************************************************************************************/
/* convert cartesian coordinates to barycentric coordinates.
   Inputs: number of cells, and the cell indices (an array),
           cartesian coordinates (an array of consecutive groups of length GEO_DIM
   Outputs: barycentric coordinates (an array of consecutive groups of length (CELL_DIM+1)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Cartesian_To_Barycentric(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             const PointType* PC,
                             PointType* PB)
{
    // first, convert from cartesian to reference domain coordinates
    PointType* PR = new PointType[Num_Cell*CELL_DIM];
    Cartesian_To_Reference(Num_Cell, CI, PC, PR);
    
    Reference_To_Barycentric(Num_Cell, PR, PB);
    delete(PR);
}

/***************************************************************************************/
/* compute the diameter of (simplex) cells in the mesh.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: diameter of each simplex (an array), defined to be the maximum edge length
   NOTE: the length of the arrays equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Diameter(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             RealType* Diam)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    if (CELL_DIM==0) // diameters are zero!
    {
        // loop through the given cells
        for (CellIndType ii=0; ii < Num_Cell; ++ii)
        {
            Diam[ii] = 0.0;
        }
    }
    else
    {
        // loop through the given cells
        for (CellIndType ii=0; ii < Num_Cell; ++ii)
        {
            // get the vertex indices of the current cell
            const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
            
            Diam[ii] = Simplex_Diameter<CELL_DIM, GEO_DIM>(c_Vtx, VI);
        }
    }
}

/***************************************************************************************/
/* return the bounding (cartesian) box of the given cells in the mesh. 
   Inputs: number of cells, and the cell indices (an array)
   Output: min and max limits of the coordinates (component-wise).
   example:  if GEO_DIM==3, then
   BB_min[] = {X_min, Y_min, Z_min},
   BB_max[] = {X_max, Y_max, Z_max}. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Bounding_Box(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             PointType* BB_min, PointType* BB_max)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const
    
    std::vector<VtxIndType> mesh_vertices;
    mesh_vertices.reserve(Num_Cell*(CELL_DIM+1));

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        for (SmallIndType tt=0; tt < (CELL_DIM+1); ++tt)
        {
            mesh_vertices.push_back(VI[tt]);
        }
    }
    
    // make the vertex list unique
    std::sort( mesh_vertices.begin(), mesh_vertices.end() );
    std::vector<VtxIndType>::iterator it;
    it = std::unique (mesh_vertices.begin(), mesh_vertices.end());
    mesh_vertices.resize( std::distance(mesh_vertices.begin(), it) );
    
    const VtxIndType Num_Vtx = (VtxIndType) mesh_vertices.size();
    const VtxIndType* VI = mesh_vertices.data();
    c_Vtx->Bounding_Box(Num_Vtx, VI, BB_min, BB_max);
}
/* return the bounding (cartesian) box of the *entire* mesh. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Bounding_Box(PointType* BB_min, PointType* BB_max)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const
    
    std::vector<VtxIndType> mesh_vertices;
    this_mesh.Get_Unique_Vertices(mesh_vertices);
    
    const VtxIndType Num_Vtx = (VtxIndType) mesh_vertices.size();
    const VtxIndType* VI = mesh_vertices.data();
    c_Vtx->Bounding_Box(Num_Vtx, VI, BB_min, BB_max);
}

/***************************************************************************************/
/* compute the volume of (simplex) cells in the mesh.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: volume of each simplex (an array)
   NOTE: the length of the arrays equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Volume(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             RealType* Vol)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        
        Vol[ii] = Simplex_Volume<CELL_DIM, GEO_DIM>(c_Vtx, VI);
    }
}

/***************************************************************************************/
/* compute the barycenter of (simplex) cells in the mesh.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: cartesian coordinates of center (an array of consecutive groups of length GEO_DIM)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Barycenter(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             PointType* BC)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);

        // get relevant output for current cell
        PointType* BC_local = BC + ii*GEO_DIM;
        
        Simplex_Barycenter<CELL_DIM, GEO_DIM>(c_Vtx, VI, BC_local);
    }
}

/***************************************************************************************/
/* compute the circumcenter and circumradius of (simplex) cells in the mesh.
   see:  https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm#Coxeter
   Inputs: number of cells, and the cell indices (an array)
   Outputs: barycentric coordinates of center (an array of consecutive groups of length (CELL_DIM+1))
            circumradius (an array with same length as number of cells)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Circumcenter(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             PointType* CB, RealType* CR)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        
        // get relevant output for current cell
        PointType* CB_ii = CB + ii*(CELL_DIM+1);
        
        CR[ii] = Simplex_Circumcenter<CELL_DIM, GEO_DIM>(c_Vtx, VI, CB_ii);
    }
}

/***************************************************************************************/
/* compute the incenter and inradius of (simplex) cells in the mesh.
   see:  "Coincidences of simplex centers and related facial structures" by Edmonds, et al.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: barycentric coordinates of center (an array of consecutive groups of length (CELL_DIM+1))
            inradius (an array with same length as number of cells)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Incenter(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             PointType* CB, RealType* CR)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);

        // get relevant output for current cell
        PointType* CB_ii = CB + ii*(CELL_DIM+1);
        
        CR[ii] = Simplex_Incenter<CELL_DIM, GEO_DIM>(c_Vtx, VI, CB_ii);
    }
}

/***************************************************************************************/
/* compute the "shape regularity" of (simplex) cells in the mesh.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: the "shape regularity" ratio (an array with same length as number of cells)
*/
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Shape_Regularity(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             RealType* SR)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);
        
        SR[ii] = Simplex_Shape_Regularity<CELL_DIM, GEO_DIM>(c_Vtx, VI);
    }
}

/***************************************************************************************/
/* compute the interior "dihedral" angles of the facets of (simplex) cells in the mesh.
   Inputs: number of cells, and the cell indices (an array)
   Outputs: the "shape regularity" ratio (an array with same length as number of cells)
   Outputs: angles in radians (an array of consecutive groups of length
              (CELL_DIM*(CELL_DIM+1)/2); this is the number of independent angles)
   NOTE: the number of groups equals the number of cells. */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MMC<CELL_DIM, GEO_DIM>::Angles(
                             const CellIndType& Num_Cell, const CellIndType* CI,
                             RealType* Ang)
{
    MMC const& this_mesh = *this; // for const
    const BasePtCoord<GEO_DIM>* const& c_Vtx = Vtx; // for const

    // loop through the given cells
    for (CellIndType ii=0; ii < Num_Cell; ++ii)
    {
        // get the vertex indices of the current cell
        const VtxIndType* VI = this_mesh.Get_Cell_vtx(CI[ii]);

        // get relevant output for current cell
        RealType* Ang_ii = Ang + ii*(CELL_DIM*(CELL_DIM+1)/2);
        
        Simplex_Angles<CELL_DIM, GEO_DIM>(c_Vtx, VI, Ang_ii);
    }
}




// SWW: the only methods that go here are the ones that require BOTH cell connectivity
//      *and* vertex coordinates.

// add these methods:

// is it better to have a generic method that calls the other stuff???

// make simple Python interface BEFORE doing mesh refinement...

/* easy:

(need a local-to-global mapping sub-routine...)

Remove_Unused_Vertices

*/

/*

function [Global_to_Local_Indexing, Local_to_Global_Indexing] =...
                 Create_Global_Local_Index_Mapping(Ordered_Indices,Max_Global_Index)
%Create_Global_Local_Index_Mapping
%
%   This routine creates two vectors of data that act as mappings between 'local' and
%   'global' indexing.
%
%   [Global_to_Local_Indexing, Local_to_Global_Indexing] =
%              Create_Global_Local_Index_Mapping(Ordered_Indices,Max_Global_Index)
%
%   OUTPUTS
%   -------
%   Global_to_Local_Indexing:
%       A (Max_Global_Index x 1) array that maps from global indices to local
%       indices.  In other words, to get the local index corresponding to global
%       index 'j', just look at Global_to_Local_Indexing(j).  Note: some of the riws
%
%   Local_to_Global_Indexing:
%       An Nx1 array that maps from local indices to global indices.  In other words,
%       to get the global index corresponding to local index 'i', just look at
%       Local_to_Global_Indexing(i).
%       Note: 'Local_to_Global_Indexing' = 'Ordered_Indices'.
%
%   Examples:
%
%   Global_to_Local_Indexing:
%
%         (Global Indexing)     (Col)       (Local Indexing)
%          Node #1 (Row #1):      6     <-- Local Node #6
%          Node #2 (Row #2):      3     <-- Local Node #3
%          Node #3 (Row #3):      0     <-- NOT a Local Node
%          Node #4 (Row #4):      0     <-- NOT a Local Node
%          Node #5 (Row #5):      2     <-- Local Node #2
%          Node #6 (Row #6):      1     <-- Local Node #1
%          Node #7 (Row #7):      0     <-- NOT a Local Node
%          Node #8 (Row #8):      7     <-- Local Node #7
%          ...
%
%   Local_to_Global_Indexing:
%
%          (Local Indexing)     (Col)      (Global Indexing)
%          Node #1 (Row #1):      6     <-- Global Node #6
%          Node #2 (Row #2):      5     <-- Global Node #5
%          Node #3 (Row #3):      2     <-- Global Node #2
%          Node #4 (Row #4):     10     <-- Global Node #10
%          Node #5 (Row #5):     14     <-- Global Node #14
%          Node #6 (Row #6):      1     <-- Global Node #1
%          Node #7 (Row #7):      8     <-- Global Node #8
%          Node #8 (Row #8):     17     <-- Global Node #17
%          ...
%
%   INPUTS
%   ------
%   Ordered_Indices:
%       An Nx1 array of numbers considered to be indices into a global data
%       structure.  However, row i of 'Ordered_Indices' corresponds to index 'i' in
%       the local data structure.
%
%   Max_Global_Index:
%       Largest index value in the global data structure (must be greater than N).

% Copyright (c) 06-05-2007,  Shawn W. Walker

Local_to_Global_Indexing = Ordered_Indices;
Local_Index_Vec = (1:1:length(Ordered_Indices))';
Global_to_Local_Indexing = zeros(Max_Global_Index,1);
Global_to_Local_Indexing(Local_to_Global_Indexing,1) = Local_Index_Vec;

end

*/

// make uniform refinement a separate thing

// need adaptive refinement

// don't worry about modifying the cell data structure yet, b/c this is connected to refinement

// need topological/connectivity changes!!!!  for this, in 3-D, can think of an edge as the intersection of two half-facets?


/* easy:

faceNormal (simplex normal/normal space)
featureEdges (complicated for general dimension)
vertexNormal (averaged normal space when normal space is 1-D)

tetrahedron:

Order_Cell_Vertices_For_Hcurl (special case for tetrahedra, CELL_DIM==3)
faces?
facetNormal?
Get_Facet_Info?

interval:

edgeTangent (simplex tangent/tangent space)
featurePoints
Get_Arclength
vertexTangent (averaged tangent space when tangent space is 1-D)

*/

/* harder:

Get_Adjacency_Matrix
Reorder
Get_Facet_Info?
Get_Laplacian_Smoothing_Matrix

**Refine

*/

/* move to multi-mesh

nearestNeighbor:  Vertex closest to a specified point
pointLocation:   simplex containing specified point

*/


// For my mesh generator, need to evaluate cut edges. can do this by first finding all intersected tetrahedra.  then make an edge mesh from that?
// Or just store all the unique edges of the mesh?  think of the edge mesh as a graph!  Store each edge as a pair of vertices (with smallest global index first) and the cell that the edge belongs to.


#undef MMC

/***/
