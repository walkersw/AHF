/*
============================================================================================
   Class for storing a vector of point coordinates.  Meant to work with BaseMesh.
   (This is just a glorofied C++ STL:vector class.)

   Note: Everything is indexed starting at 0!

   Copyright (c) 05-18-2020,  Shawn W. Walker
============================================================================================
*/

#define _BASEPTCOORD_CC

#ifndef _PRELIM_H
#include "Prelim.h" // basic typedefs and includes
#endif

#ifndef _BASICSTRUCTS_H
#include "BasicStructs.h" // basic structs for the BasePtCoord class
#endif

/* C++ class definition */
#define  BPC  BasePtCoord
// template the geometric dimension (of the vertex coordinates)
template <SmallIndType GEO_DIM>
class BPC
{
public:
    BPC();
    ~BPC();

    // open Point for modification
    inline void Open() { Point_Open = true; };
    // close Point; modification is no longer allowed
    inline void Close() { Point_Open = false; };
    // generic check for if mesh is open for modification
    inline bool Is_Point_Open() const
    {
        if (!Point_Open)
        {
            std::cout << "Point is not open for modification!" << std::endl;
            std::cout << "     You must first use the 'Open' method." << std::endl;
        }
        return Point_Open;
    };

    void Clear() // clear all data
    {
        Point.clear();
        Point_Open = true; // re-open Point for modification
    };

    // get the topological dimension
    inline SmallIndType Geo_Dim() const { return GEO_DIM; };
    // get number of vertex coordinates stored
    inline VtxIndType Num_Points() const { return (VtxIndType) Point.size(); };
    // allocate room for specified number of vertex point coordinates (plus some room to grow)
    void Reserve_Points(const VtxIndType&);
    // initialize point coordinates to the *origin* for a specific number of vertices
    void Init_Points(const VtxIndType&);
    
    // set all of the point coordinates at once
    void Set_Coord_Data(const PointType*, const VtxIndType&);
    // set a subset of point coordinates
    void Set_Coord_Data(const PointType*, const VtxIndType&, const VtxIndType&);
    // set point coordinates of specific vertex
    void Set_Coord(const VtxIndType&, const PointType*);
    void Set_Coord(const VtxIndType&, const PointType&);
    void Set_Coord(const VtxIndType&, const PointType&, const PointType&);
    void Set_Coord(const VtxIndType&, const PointType&, const PointType&, const PointType&);

    // retrieve one point (writeable)
    inline PointType* Get_Point_coord(const VtxIndType&);
    inline VtxCoordType<GEO_DIM>& Get_Point(const VtxIndType&);

    // retrieve one point (read-only)
    inline const PointType* Get_Point_coord(const VtxIndType&) const;
    inline const VtxCoordType<GEO_DIM>& Get_Point(const VtxIndType&) const;

    // get the cartesian box that bounds all the points
    void Bounding_Box(PointType*, PointType*) const;
    void Bounding_Box(const VtxIndType&, const VtxIndType*, PointType*, PointType*) const;
    
    // print out vertex coordinates
    void Display_Vtx_Coord(const VtxIndType& vi=NULL_Vtx) const;

protected:
    /* main data storage */
    typedef VtxCoordType<GEO_DIM> VtxCoord_DIM; // convenient
    std::vector<VtxCoord_DIM>  Point; // store vertex point coordinates here

    double   Point_Reserve_Buffer; // amount of extra memory to allocate when re-allocating
                                   // Point

    /* flag to indicate if points may be added or changed.
       true  = point coordinates may be changed, etc.
       false = the points cannot be changed! */
    bool Point_Open;
};

/***************************************************************************************/
/* constructor */
template <SmallIndType GEO_DIM>
BPC<GEO_DIM>::BPC()
{
    Point_Open = true; // Point starts out as open for modification

    // ensure memory is clear to start
    Clear();
    Point_Reserve_Buffer = 0.2; // allocate an extra 20% when re-allocating vertex coordinates

    //std::cout << "BasePtCoord constructor..." << std::endl;
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType GEO_DIM>
BPC<GEO_DIM>::~BPC()
{
    // clear the data
    Clear();

    //std::cout << "BasePtCoord destructor..." << std::endl;
}

/***************************************************************************************/
/* Allocate memory to hold vertex point coordinate data (plus a little). */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Reserve_Points(const VtxIndType& Num_Pts)
{
    if (!Is_Point_Open())
        return;

    // compute the actual size to allocate for the points
    const VtxIndType Actual_Point_SIZE = (VtxIndType) ((1.0 + Point_Reserve_Buffer) * Num_Pts);
    Point.reserve(Actual_Point_SIZE);
}

/***************************************************************************************/
/* Initialize Point to contain Num_Pts vertices, whose coordinates are initialized
   to the origin (0,0,...,0). */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Init_Points(const VtxIndType& Num_Pts)
{
    if (!Is_Point_Open())
        return;

    VtxCoord_DIM Zero_Pt;
    Zero_Pt.Set(); // init to (0,0,...,0)

    if (Point.capacity() < Num_Pts)
    {
        Reserve_Points(Num_Pts);
    }
    else
    {
        Point.resize(Num_Pts);
    }
	
	typename std::vector<VtxCoord_DIM>::iterator it;
    it = Point.begin();
    Point.insert(it,Num_Pts,Zero_Pt); // fill in
}

/***************************************************************************************/
/* Set all the coordinates at once.
   Note: input is the global coordinates of all the vertices. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord_Data(const PointType* Coord_Data, const VtxIndType& Num_Coord_Data)
{
    Init_Points(Num_Coord_Data);
    
    // now fill it in
    for (VtxIndType ii = 0; ii < Num_Coord_Data; ++ii)
    {
        Point[ii].Set(Coord_Data + GEO_DIM*ii);
    }
}

/***************************************************************************************/
/* Set a subset of the coordinates.
   Note: input is the global coordinates of subset of vertices. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord_Data(const PointType* Coord_Data,
                                  const VtxIndType& Begin_Index, const VtxIndType& Num_Pts)
{
    if (!Is_Point_Open())
        return;

    const VtxIndType Actual_Num_Pts = Begin_Index + Num_Pts;
    if (Point.capacity() < Actual_Num_Pts)
    {
        Reserve_Points(Actual_Num_Pts);
    }
    if (Point.size() < Actual_Num_Pts)
    {
        VtxCoord_DIM Zero_Pt;
        Zero_Pt.Set(); // init to (0,0,...,0)
        Point.resize(Actual_Num_Pts,Zero_Pt);
    }
    
    // now fill it in
    for (VtxIndType ii = Begin_Index; ii < Actual_Num_Pts; ++ii)
    {
        Point[ii].Set(Coord_Data + GEO_DIM*ii);
    }
}

/***************************************************************************************/
/* Set the coordinates of a specific vertex. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi, const PointType* vtx_coord)
{
    if (!Is_Point_Open())
        return;

    assert( vi < Num_Points() );
    Point[vi].Set(vtx_coord);
}
// 1-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0)
{
    if (!Is_Point_Open())
        return;

    assert( vi < Num_Points() );
    Point[vi].Set(x0);
}
// 2-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0, const PointType& x1)
{
    if (!Is_Point_Open())
        return;

    assert( vi < Num_Points() );
    Point[vi].Set(x0,x1);
}
// 3-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0, const PointType& x1, const PointType& x2)
{
    if (!Is_Point_Open())
        return;

    assert( vi < Num_Points() );
    Point[vi].Set(x0,x1,x2);
}

/***************************************************************************************/
/* Get reference to specific mesh vertex coordinate data (given the vertex index). */
// WRITEABLE
template <SmallIndType GEO_DIM>
inline PointType* BPC<GEO_DIM>::Get_Point_coord(const VtxIndType& vi) // writeable
{
    if (!Is_Point_Open())
    {
        std::cerr << "Fatal error in 'Get_Point_coord'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // vi must be in [0, Num_Points), and not invalid
    assert((vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi].coord;
}
/* Get pointer to specific point data (given the point vertex index). */
template <SmallIndType GEO_DIM>
inline VtxCoordType<GEO_DIM>& BPC<GEO_DIM>::Get_Point(const VtxIndType& vi) // writeable
{
    if (!Is_Point_Open())
    {
        std::cerr << "Fatal error in 'Get_Point'!" << std::endl;
        std::cerr << "     Mesh is not 'open' for writing." << std::endl;
        std::exit(1);
    }

    // vi must be in [0, Num_Points), and not invalid
    assert((vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi];
}
// READ-ONLY
template <SmallIndType GEO_DIM>
inline const PointType* BPC<GEO_DIM>::Get_Point_coord(const VtxIndType& vi) const // read-only
{
    // vi must be in [0, Num_Points), and not invalid
    assert((vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi].coord;
}
/* Get pointer to specific point data (given the point vertex index). */
template <SmallIndType GEO_DIM>
inline const VtxCoordType<GEO_DIM>& BPC<GEO_DIM>::Get_Point(const VtxIndType& vi) const // read-only
{
    // vi must be in [0, Num_Points), and not invalid
    assert((vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi];
}

/***************************************************************************************/
/* Get the bounding box of all the vertex coordinates.
   Output: min and max limits of the coordinates (component-wise).
   example:  if GEO_DIM==3, then
   BB_min[] = {X_min, Y_min, Z_min},
   BB_max[] = {X_max, Y_max, Z_max}. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Bounding_Box(PointType* BB_min, PointType* BB_max) const
{
    // initialize the box with the 0th point
    const VtxCoord_DIM& P0 = Point[0];
    for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
    {
        BB_min[ii] = P0.coord[ii];
        BB_max[ii] = P0.coord[ii];
    }
    
    // now loop through the rest...
    for (std::vector<VtxCoord_DIM>::iterator it = Point.begin(); it != Point.end(); ++it)
    {
        const VtxCoord_DIM& PT = *it;
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
        {
            if (PT.coord[ii] < BB_min[ii]) BB_min[ii] = PT.coord[ii];
            if (PT.coord[ii] > BB_max[ii]) BB_max[ii] = PT.coord[ii];
        }
    }
}
/* this allows to only look through a subset of vertex indices in "VI". */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Bounding_Box(const VtxIndType& Num_Vtx, const VtxIndType* VI,
                                PointType* BB_min, PointType* BB_max) const
{
    // initialize the box with the 0th point
    const VtxCoord_DIM& P0 = Point[VI[0]];
    for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
    {
        BB_min[ii] = P0.coord[ii];
        BB_max[ii] = P0.coord[ii];
    }
    
    // now loop through the rest...
    for (VtxIndType vi = 0; vi < Num_Vtx; ++vi)
    {
        const VtxCoord_DIM& PT = Point[VI[vi]];
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
        {
            if (PT.coord[ii] < BB_min[ii]) BB_min[ii] = PT.coord[ii];
            if (PT.coord[ii] > BB_max[ii]) BB_max[ii] = PT.coord[ii];
        }
    }
}

/***************************************************************************************/
/* print vertex coordinates. "vi" is the index of a specific vertex;
   if vi=NULL_Vtx, then print all vertex coordinates. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Display_Vtx_Coord(const VtxIndType& main_vi) const
{
    if (main_vi==NULL_Vtx)
    {
        // then print all the vertices' coordinates
        std::cout << "Display coordinates of all vertices:" << std::endl;
        std::cout << "Vtx #    |        Coordinates" << std::endl;
        for (VtxIndType vi = 0; vi < Num_Points(); ++vi)
        {
            const PointType* VC = Get_Point_coord(vi);
            std::cout << vi << "  |  " << "(" << VC[0]; // print vtx # and first coordinate value
            for (SmallIndType kk = 1; kk < GEO_DIM; ++kk)
                std::cout << ", " << VC[kk]; // print other coordinates
            std::cout << ")" << std::endl;
        }
    }
    else
    {
        // then print ONE vertex
        std::cout << "Display coordinates of vertex #" << main_vi << ":" << std::endl;

        const PointType* VC = Get_Point_coord(main_vi);
        std::cout << "  " << "(" << VC[0]; // print vtx # and first coordinate value
        for (SmallIndType kk = 1; kk < GEO_DIM; ++kk)
            std::cout << ", " << VC[kk]; // print other coordinates
        std::cout << ")" << std::endl;
    }
}

// SWW: what else do we need here?

#undef BPC

/***/
