/*
============================================================================================
   Definition of small light-weight structs used to hold single mesh cell data, single
   vertex coordinate data, etc.

   Note: these structs are used within BaseMesh and BasePtCoord, in the form of
         std::vectors of cells and point coordinates.
   Note: everything is indexed from 0!

   Copyright (c) 12-15-2016,  Shawn W. Walker
============================================================================================
*/

#define _BASICSTRUCTS_H

#ifndef _PRELIM_H
#include "Prelim.h" // basic typedefs and includes
#endif
#ifndef _VTX2HALFFACET_MAPPING_CC
#include "Vtx2HalfFacet_Mapping.cc" // class for handling the Vertex-to-Half-Facet Mapping
#endif

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
struct CellType // this is an interface
{
    // get internal variables
    virtual VtxIndType* Get_vtx()=0;
    virtual HalfFacetType* Get_halffacet()=0;
    virtual const VtxIndType* Get_vtx() const = 0;
    virtual const HalfFacetType* Get_halffacet() const = 0;

    // check equality
    //virtual inline bool Equal(const CellType*) const = 0;
    virtual bool Equal(const CellType&) const = 0;
    virtual bool Equal(const VtxIndType*, const HalfFacetType*) const = 0;

    // set *one* vertex/facet index
    virtual CellType& Set(const SmallIndType&, const VtxIndType&, const HalfFacetType&)=0;
    // set all vertex indices, and init half-facets to NULL
    virtual CellType& Set(const VtxIndType*)=0;
    // 1-D
    virtual CellType& Set(const VtxIndType&, const VtxIndType&)=0;
    // 2-D
    virtual CellType& Set(const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    // 3-D
    virtual CellType& Set(const VtxIndType&, const VtxIndType&, const VtxIndType&, const VtxIndType&)=0;
    // simple print command
    virtual void Print() const = 0;
};
/*** here is the actual concrete base class ***/
template<SmallIndType CELL_DIM>
struct CellSimplexType : CellType
{
    VtxIndType           vtx[CELL_DIM+1]; // global vertex indices of the cell
    HalfFacetType  halffacet[CELL_DIM+1]; // half-facets corresponding to local facets of cell
                                          // note: see examples in concrete sub-class
    // get internal variables
    inline VtxIndType* Get_vtx() { return vtx; };
    inline HalfFacetType* Get_halffacet() { return halffacet; };
    inline const VtxIndType* Get_vtx() const { return vtx; };
    inline const HalfFacetType* Get_halffacet() const { return halffacet; };

    // check equality
    // inline bool Equal(const CellType* IN) const
    // {
        // return Equal(*IN);
    // }
    inline bool Equal(const CellType& IN) const
    {
        const VtxIndType* IN_vtx = IN.Get_vtx();
        const HalfFacetType* IN_halffacet = IN.Get_halffacet();
        for (SmallIndType ii = 0; ii < (CELL_DIM+1); ++ii)
        {
            if (IN_vtx[ii]!=vtx[ii]) return false;
            if (!IN_halffacet[ii].Equal(halffacet[ii])) return false;
        }
        return true;
    }
    inline bool Equal(const CellSimplexType<CELL_DIM>& IN) const
    {
        for (SmallIndType ii = 0; ii < (CELL_DIM+1); ++ii)
        {
            if (IN.vtx[ii]!=vtx[ii]) return false;
            if (!IN.halffacet[ii].Equal(halffacet[ii])) return false;
        }
        return true;
    }
    inline bool Equal(const VtxIndType* IN_vtx, const HalfFacetType* IN_halffacet) const
    {
        for (SmallIndType ii = 0; ii < (CELL_DIM+1); ++ii)
        {
            if (IN_vtx[ii]!=vtx[ii]) return false;
            if (!IN_halffacet[ii].Equal(halffacet[ii])) return false;
        }
        return true;
    }
    // set *one* vertex/facet index
    inline CellType& Set(const SmallIndType& index, const VtxIndType& v, const HalfFacetType& hf)
    {
        vtx[index] = v;
        if ( ( (hf.fi > CELL_DIM) ) && !hf.Is_Null() )
        {
            std::cout << "ERROR: facet index fi is too small or too large!" << std::endl;
            std::cout << "       Setting half-facet to NULL..." << std::endl;
            halffacet[index].Set();
        }
        else
            halffacet[index].Set(hf);

        return *this;
    }
    // set all vertex indices, and init half-facets to NULL
    inline CellType& Set(const VtxIndType* v)
    {
        for (SmallIndType ii = 0; ii < (CELL_DIM+1); ++ii)
        {
            vtx[ii] = v[ii];
            // init to null value
            halffacet[ii].Set();
        }
        return *this;
    }
    // 1-D
    inline CellType& Set(const VtxIndType& v0, const VtxIndType& v1)
    {
        assert(CELL_DIM==1);
        vtx[0] = v0;
        vtx[1] = v1;
        halffacet[0].Set();
        halffacet[1].Set();

        return *this;
    }
    // 2-D
    inline CellType& Set(const VtxIndType& v0, const VtxIndType& v1, const VtxIndType& v2)
    {
        assert(CELL_DIM==2);
        vtx[0] = v0;
        vtx[1] = v1;
        vtx[2] = v2;
        halffacet[0].Set();
        halffacet[1].Set();
        halffacet[2].Set();

        return *this;
    }
    // 3-D
    inline CellType& Set(const VtxIndType& v0, const VtxIndType& v1,
                         const VtxIndType& v2, const VtxIndType& v3)
    {
        assert(CELL_DIM==3);
        vtx[0] = v0;
        vtx[1] = v1;
        vtx[2] = v2;
        vtx[3] = v3;
        halffacet[0].Set();
        halffacet[1].Set();
        halffacet[2].Set();
        halffacet[3].Set();

        return *this;
    }
    // simple print command
    inline void Print() const
    {
        std::cout << "Cell Vertex Indices: " << vtx[0];
        for (SmallIndType ii = 1; ii < (CELL_DIM+1); ++ii)
            std::cout << ", " << vtx[ii];
        std::cout << std::endl;
        std::cout << "Cell Half-facets: ";
        halffacet[0].Print();
        for (SmallIndType ii = 1; ii < (CELL_DIM+1); ++ii)
        {
            std::cout << ", ";
            halffacet[ii].Print();
        }
        std::cout << std::endl;
    };
};

/***************************************************************************************/
/* REMARK: We need to define a one-to-one mapping between *local* facets and
   *local* vertices of a cell.

   For an edge (CELL_DIM==1), we define this to be:
      Vtx | Facet (Vertex)
     -----+-------------
       0  |   0
       1  |   1

   For a triangle (CELL_DIM==2), we define this to be:
      Vtx | Facet (Edge)
     -----+-------------
       0  |   1
       1  |   2
       2  |   0

   For a tetrahedron (CELL_DIM==3), we define this to be:
      Vtx | Facet (Face)
     -----+-------------
       0  |   1
       1  |   2
       2  |   3
       3  |   0

   (The pattern is obvious for higher cell dimensions...)

   Note: each vertex is contained (attached) to its corresponding facet.
*/

/*****************************************************************************/
/* struct for "marking" which vertices/half-facets have been visited. */
template<SmallIndType CELL_DIM>
struct CellSimplexMarkedType
{
    bool      facet[CELL_DIM+1]; // local facets (true/false = visited/or not)
};

/* struct for holding a single edge of a mesh.
   An edge is defined by [v0, v1], where v0, v1 are the end vertices of the edge.
   Note: for a simplex mesh, the edge [v0, v1] exists when v0, v1 are both
         contained in the *same* mesh cell. */
struct MeshEdgeType
{
    VtxIndType      vtx[2]; // global vertex indices of end points of an edge
	
    // check equality
    inline bool Equal(const MeshEdgeType& IN) const
    {
        return ( (IN.vtx[0]==vtx[0]) && (IN.vtx[1]==vtx[1]) );
    }
    // check for NULL
    inline bool Is_Null() const
    {
        if (vtx[0]==NULL_Vtx || vtx[1]==NULL_Vtx)
            return true;
        else
            return false;
    }
    // set the two vertices
    inline void Set(const VtxIndType& v0, const VtxIndType& v1)
    {
        vtx[0] = v0;
        vtx[1] = v1;
    }
    // sort the two vertices (ascending order)
    inline void Sort()
    {
        const VtxIndType& temp_v0 = vtx[1];
        const VtxIndType& temp_v1 = vtx[0];
        if (vtx[0] > vtx[1])
            Set(temp_v0, temp_v1);
        // otherwise, already sorted
    }
    // set the vertices in sorted (ascending) order
    inline void Set_Sorted(const VtxIndType& v0, const VtxIndType& v1)
    {
        if (v0 < v1)
            Set(v0, v1);
        else
            Set(v1, v0);
    }
    // simple print command
    inline void Print() const
    {
        if (Is_Null())
            std::cout << "[" << "-" << ", " << "-" << "]";
        else
            std::cout << "[" << vtx[0] << ", " << vtx[1] << "]";
    };
};

/***************************************************************************************/
/* vertex coordinate data.
   GEO_DIM is the geometric dimension of the coordinates.
*/
struct CoordType // this is an interface
{
    // get internal variables
    virtual PointType* Get_coord()=0;
    virtual const PointType* Get_coord() const = 0;

    // check equality
    virtual bool Equal(const CoordType&) const = 0;
    // set point coordinates
    virtual CoordType& Set(const PointType*)=0;
    // 1-D
    virtual CoordType& Set(const PointType&)=0;
    // 2-D
    virtual CoordType& Set(const PointType&, const PointType&)=0;
    // 3-D
    virtual CoordType& Set(const PointType&, const PointType&, const PointType&)=0;
    // default: origin
    virtual CoordType& Set()=0;
    // print
    virtual void Print() const = 0;
};
/*** here is the actual concrete base class ***/
template<SmallIndType GEO_DIM>
struct VtxCoordType : CoordType
{
    PointType     coord[GEO_DIM]; // vertex coordinates

    // get internal variables
    inline PointType* Get_coord() { return coord; };
    inline const PointType* Get_coord() const { return coord; };

    // check equality
    inline bool Equal(const CoordType& IN) const
    {
        const PointType* IN_coord = IN.Get_coord();
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            if (IN_coord[ii]!=coord[ii]) return false;

        return true;
    }
    inline bool Equal(const VtxCoordType<GEO_DIM>& IN) const
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            if (IN.coord[ii]!=coord[ii]) return false;

        return true;
    }
    // set point coordinates
    inline CoordType& Set(const PointType* p)
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            coord[ii] = p[ii];

        return *this;
    }
    // 1-D
    inline CoordType& Set(const PointType& x0)
    {
        assert(GEO_DIM==1);
        coord[0] = x0;
        return *this;
    }
    // 2-D
    inline CoordType& Set(const PointType& x0, const PointType& x1)
    {
        assert(GEO_DIM==2);
        coord[0] = x0;
        coord[1] = x1;
        return *this;
    }
    // 3-D
    inline CoordType& Set(const PointType& x0, const PointType& x1, const PointType& x2)
    {
        assert(GEO_DIM==3);
        coord[0] = x0;
        coord[1] = x1;
        coord[2] = x2;
        return *this;
    }
    // default: origin
    inline CoordType& Set()
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            coord[ii] = 0.0;

        return *this;
    }
    // simple print command
    inline void Print() const
    {
        std::cout << "Vertex coordinates: (" << coord[0];
        for (SmallIndType ii = 1; ii < GEO_DIM; ++ii)
            std::cout << ", " << coord[ii];
        std::cout << ")" << std::endl;
    };
};

// SWW: can maybe have some more convenience routines here...

/***/
