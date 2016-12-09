/*
============================================================================================
   Class for storing a vector of point coordinates.

   Note: Everything is indexed starting at 0!

   Copyright (c) 12-08-2016,  Shawn W. Walker
============================================================================================
*/

#define _BASEPTCOORD_CC

#ifndef _PRELIM_H
#include "Prelim.h" // basic typedefs and includes
#endif

/***************************************************************************************/
/* vertex coordinate data.
   GEO_DIM is the geometric dimension of the coordinates.
*/
template<SmallIndType GEO_DIM>
struct VtxCoordType
{
    PointType     coord[GEO_DIM]; // vertex coordinates
	
    // check equality
    inline bool Equal(const VtxCoordType<GEO_DIM>& IN) const
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            if (IN.coord[ii]!=coord[ii]) return false;

        return true;
    }
    // set point coordinates
    inline VtxCoordType<GEO_DIM>& Set(const PointType p[GEO_DIM])
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            coord[ii] = p[ii];

		return *this;
    }
	// 1-D
	inline VtxCoordType<GEO_DIM>& Set(const PointType x0)
    {
		assert(GEO_DIM==1);
        coord[0] = x0;
		return *this;
    }
	// 2-D
	inline VtxCoordType<GEO_DIM>& Set(const PointType x0, const PointType x1)
    {
		assert(GEO_DIM==2);
        coord[0] = x0;
		coord[1] = x1;
		return *this;
    }
	// 3-D
	inline VtxCoordType<GEO_DIM>& Set(const PointType x0, const PointType x1, const PointType x2)
    {
		assert(GEO_DIM==3);
        coord[0] = x0;
		coord[1] = x1;
		coord[2] = x2;
		return *this;
    }
    // default: origin
    inline VtxCoordType<GEO_DIM>& Set()
    {
        for (SmallIndType ii = 0; ii < GEO_DIM; ++ii)
            coord[ii] = 0.0;

		return *this;
    }
};

/* C++ class definition */
#define  BPC  BasePtCoord
// template the geometric dimension (of the vertex coordinates)
template <SmallIndType GEO_DIM>
class BPC
{
public:
    BPC();
    ~BPC();
    void Clear() // clear all data
    {
		Point.clear();
    };
	
	typedef VtxCoordType<GEO_DIM> VtxCoord_DIM; // convenient
	
	// get the topological dimension
	inline SmallIndType Geo_Dim() const { return GEO_DIM; };
    // get number of vertex coordinates stored
    inline VtxIndType Num_Points() const { return (VtxIndType) Point.size(); };
	// allocate room for specified number of vertex point coordinates (plus some room to grow)
    void Reserve_Points(const VtxIndType&);
	// initialize point coordinates to the *origin* for a specific number of vertices
    void Init_Points(const VtxIndType&);
	// set point coordinates of specific vertex
    void Set_Coord(const VtxIndType& vi, const PointType vc[GEO_DIM]);
	void Set_Coord(const VtxIndType&, const PointType&);
	void Set_Coord(const VtxIndType&, const PointType&, const PointType&);
	void Set_Coord(const VtxIndType&, const PointType&, const PointType&, const PointType&);
	// retrieve one point
    inline VtxCoord_DIM& Get_Point(const VtxIndType&); // for writeable
    inline const VtxCoord_DIM& Get_Point(const VtxIndType&) const; // for read-only

    // print out vertex coordinates
    void Display_Vtx_Coord(const VtxIndType& vi=NULL_Vtx) const;
	
	/* main data storage */
	std::vector<VtxCoord_DIM>  Point; // store vertex point coordinates here

protected:
    double   Point_Reserve_Buffer; // amount of extra memory to allocate when re-allocating
                                   // Point
};

/***************************************************************************************/
/* constructor */
template <SmallIndType GEO_DIM>
BPC<GEO_DIM>::BPC()
{
    // ensure memory is clear to start
    Clear();
	Point_Reserve_Buffer = 0.2; // allocate an extra 20% when re-allocating vertex coordinates
}

/***************************************************************************************/
/* DE-structor */
template <SmallIndType GEO_DIM>
BPC<GEO_DIM>::~BPC()
{
    // clear the data
    Clear();
}

/***************************************************************************************/
/* Allocate memory to hold vertex point coordinate data (plus a little). */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Reserve_Points(const VtxIndType& Num_Pts)
{
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
	VtxCoordType<GEO_DIM> Zero_Pt;
	Zero_Pt.Set(); // init to (0,0,...,0)
	
	if (Point.capacity() < Num_Pts)
		Reserve_Points(Num_Pts);

	std::vector<VtxCoord_DIM>::iterator it;
	it = Point.begin();
	Point.insert(it,Num_Pts,Zero_Pt); // fill in
}

/***************************************************************************************/
/* Set the coordinates of a specific vertex. */
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi, const PointType vtx_coord[GEO_DIM])
{
	assert( vi < Num_Points() );
	Point[vi].Set(vtx_coord);
}
// 1-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0)
{
	assert( vi < Num_Points() );
	Point[vi].Set(x0);
}
// 2-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0, const PointType& x1)
{
	assert( vi < Num_Points() );
	Point[vi].Set(x0,x1);
}
// 3-D
template <SmallIndType GEO_DIM>
void BPC<GEO_DIM>::Set_Coord(const VtxIndType& vi,
                             const PointType& x0, const PointType& x1, const PointType& x2)
{
	assert( vi < Num_Points() );
	Point[vi].Set(x0,x1,x2);
}

/***************************************************************************************/
/* Get reference to specific mesh vertex coordinate data (given the vertex index). */
template <SmallIndType GEO_DIM>
inline VtxCoordType<GEO_DIM>& BPC<GEO_DIM>::Get_Point(const VtxIndType& vi)
{
    // vi must be in [0, Num_Points), and not invalid
    assert((vi >= 0) && (vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi];
}
template <SmallIndType GEO_DIM>
inline const VtxCoordType<GEO_DIM>& BPC<GEO_DIM>::Get_Point(const VtxIndType& vi) const // read-only
{
    // vi must be in [0, Num_Points), and not invalid
    assert((vi >= 0) && (vi < Num_Points()) && (vi!=NULL_Vtx));

    return Point[vi];
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
            const VtxCoord_DIM& VTX = Get_Point(vi);
            std::cout << vi << "  |  " << "(" << VTX.coord[0]; // print vtx # and first coordinate value
            for (SmallIndType kk = 1; kk < GEO_DIM; ++kk)
                std::cout << ", " << VTX.coord[kk]; // print other coordinates
            std::cout << ")" << std::endl;
        }
    }
    else
    {
        // then print ONE vertex
        std::cout << "Display coordinates of vertex #" << main_vi << ":" << std::endl;

		const VtxCoord_DIM& VTX = Get_Point(main_vi);
		std::cout << "  " << "(" << VTX.coord[0]; // print vtx # and first coordinate value
		for (SmallIndType kk = 1; kk < GEO_DIM; ++kk)
			std::cout << ", " << VTX.coord[kk]; // print other coordinates
		std::cout << ")" << std::endl;
    }
}





#undef BPC

/***/
