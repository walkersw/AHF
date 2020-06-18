/*
============================================================================================
   Some simple mathematical sub-routines.

   Copyright (c) 05-23-2020,  Shawn W. Walker
============================================================================================
*/

// std C++ classes
//#include <algorithm>
//#include <vector>
//#include <map>
#include <numeric>

#include <Eigen/Dense>

#define _BASICMATH_H

#ifndef _PRELIM_H
#include "Prelim.h"  // simple typedefs, etc...
#endif
#ifndef _BASEPTCOORD_CC
#include "BasePtCoord.cc"  // base class for all vertex coordinates
#endif

static const double         PI = 3.14159265358979323846264338327950288419716939;
static const double     SQRT_2 = 1.4142135623730950488016887242097;
static const double INV_SQRT_2 = 0.7071067811865475244008443621048;
static const double     SQRT_3 = 1.7320508075688772935274463415059;
static const double INV_SQRT_3 = 0.57735026918962576450914878050196;

/* define simple routines */

/***************************************************************************************/
/* absolute value */
inline RealType my_abs(const RealType& value)
{
    if (value < 0)
        return -value;
    else
        return value;
}

/***************************************************************************************/
/* factorial */
inline unsigned int my_factorial(unsigned int n)
{
    return ((n==0) || (n==1)) ? 1 : n * my_factorial(n-1);
}

/* define mathematical sub-routines, templated with the topological dimension TD,
   and/or the geometric dimension GD. */

/***************************************************************************************/
/* get the (jacobian) matrix and translation vector for the affine map from the
   "standard" reference simplex. This is the usual finite element affine map. */
template <SmallIndType TD, SmallIndType GD>
inline void Simplex_Affine_Map(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                               Eigen::Matrix<PointType, GD, TD>& A,
                               Eigen::Matrix<PointType, GD, 1>&  b)
{
    // translation vector
    const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
    for (SmallIndType kk=0; kk < GD; ++kk)
    {
        b(kk) = XC_0[kk];
    }
    // Jacobian matrix
    for (SmallIndType tt=1; tt < (TD+1); ++tt)
    {
        // get one of the point coordinates
        const PointType* XC_tt = VX->Get_Point_coord(VI[tt]);
        for (SmallIndType gg=0; gg < GD; ++gg)
        {
            A(gg,tt-1) = XC_tt[gg] - b(gg);
        }
    }
}

/***************************************************************************************/
/* convert reference element coordinates to cartesian coordinates, where the reference
   element is the "standard" reference simplex.
   Inputs: vertices of the simplex, reference coordinates
   Outputs: cartesian coordinates */
template <SmallIndType TD, SmallIndType GD>
void Simplex_Reference_To_Cartesian(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                    const PointType* PR,
                                    PointType* PC)
{
    // compute 0th barycentric coordinate
    const PointType SUM = std::accumulate(PR,PR+TD,0.0);
    const PointType BC0 = 1.0 - SUM;

    /* note: this is computed using (effectively) barycentric coordinates */

    // init to 0th contribution
    const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
    for (SmallIndType qq=0; qq < GD; ++qq)
        PC[qq] = BC0 * XC_0[gg];

    for (SmallIndType tt=1; tt < (TD+1); ++tt)
    {
        // get one of the point coordinates
        const PointType* XC_tt = VX->Get_Point_coord(VI[tt]);
        for (SmallIndType gg=0; gg < GD; ++gg)
        {
            // compute contribution (from other barycentric coordinates)
            PC[gg] += PR[tt-1] * XC_tt[gg];
        }
    }
}

/***************************************************************************************/
/* convert cartesian coordinates to reference element coordinates, where the reference
   element is the "standard" reference simplex.
   Inputs: vertices of the simplex, cartesian coordinates
   Outputs: reference coordinates */
template <SmallIndType TD, SmallIndType GD>
void Simplex_Cartesian_To_Reference(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                    const PointType* PC,
                                    PointType* PR)
{
    // get the relevant cart-coord for the current cell
    Eigen::Matrix<PointType, GD, 1>  xc;
    for (SmallIndType gg=0; gg < GD; ++gg)
    {
        xc(gg) = PC[gg];
    }

    /* setup the affine map parameters */
    Eigen::Matrix<PointType, GD, TD>  A;
    Eigen::Matrix<PointType, GD, 1>  b;
    Simplex_Affine_Map<TD,GD>(VX, VI, A, b);

    // /* SWW: for testing: */
    // Eigen::Matrix<PointType, GD, TD> TS;
    // const RealType A_rank = Simplex_Tangent_Space<TD, GD>(VX, VI, TS);
    // // Eigen::Matrix<PointType, GD, (GD - TD)> NS;
    // // const RealType A_rank = Simplex_Normal_Space<TD, GD>(VX, VI, NS);

    // compute A' * (xc - b)
    Eigen::Matrix<PointType, TD, 1> r = A.transpose() * (xc - b);
    // compute M = A'*A
    Eigen::Matrix<PointType, TD, TD> M = A.transpose() * A;
    // solve for reference coordinates
    Eigen::Matrix<PointType, TD, 1> out = M.llt().solve(r);

    // copy into output array
    for (SmallIndType qq=0; qq < TD; ++qq)
    {
        PR[qq] = out(qq);
    }
}

/***************************************************************************************/
/* convert barycentric coordinates to reference element coordinates.
   Inputs: barycentric coordinates
   Outputs: reference coordinates */
template <SmallIndType TD>
void Simplex_Barycentric_To_Reference(const PointType* PB,
                                      PointType* PR)
{
    // go through each component
    for (SmallIndType tt=0; tt < TD; ++tt)
    {
        // note: t-th reference coordinate is the (t+1)st barycentric coordinate
        PR[tt] = PB[tt+1];
    }
}

/***************************************************************************************/
/* convert reference element coordinates to barycentric coordinates.
   Inputs: reference coordinates
   Outputs: barycentric coordinates */
template <SmallIndType TD>
void Simplex_Reference_To_Barycentric(const PointType* PR,
                                      PointType* PB)
{
    // go through each component
    PointType SUM = 0.0;
    for (SmallIndType tt=0; tt < TD; ++tt)
    {
        const PointType v = PR[tt];
        SUM += v;
        // note: t-th reference coordinate is the (t+1)st barycentric coordinate
        PB[tt+1] = v;
    }
    // the 0th barycentric coordinate is 1 - SUM
    PB[0] = 1.0 - SUM;
}

/***************************************************************************************/
/* convert barycentric coordinates to cartesian coordinates.
   Inputs: vertices of the simplex, barycentric coordinates
   Outputs: cartesian coordinates */
template <SmallIndType TD, SmallIndType GD>
void Simplex_Barycentric_To_Cartesian(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                      const PointType* PB,
                                      PointType* PC)
{
    // init to 0th contribution
    const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
    for (SmallIndType qq=0; qq < GD; ++qq)
        PC[qq] = PB[0] * XC_0[qq];

    // get the rest
    for (SmallIndType tt=1; tt < (TD+1); ++tt)
    {
        // get one of the other point coordinates
        const PointType* XC_tt = VX->Get_Point_coord(VI[tt]);
        for (SmallIndType gg=0; gg < GD; ++gg)
        {
            // compute contribution
            PC[gg] += PB[tt] * XC_tt[gg];
        }
    }
}

/***************************************************************************************/
/* convert cartesian coordinates to barycentric coordinates.
   Inputs: vertices of the simplex, cartesian coordinates
   Outputs: barycentric coordinates */
template <SmallIndType TD, SmallIndType GD>
void Simplex_Cartesian_To_Barycentric(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                      const PointType* PC,
                                      PointType* PB)
{
    // first, convert from cartesian to reference domain coordinates
    PointType PR[TD];
    Simplex_Cartesian_To_Reference<TD,GD>(VX, VI, PC, PR);
    Simplex_Reference_To_Barycentric<TD>(PR, PB);
}

/***************************************************************************************/
/* compute the diameter of a simplex */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Diameter(const BasePtCoord<GD>* const& VX, const VtxIndType* VI)
{
    // compute all edge lengths, and take the MAX
    RealType MAX_Edge_Length = 0.0;
    for (SmallIndType rr=0; rr < (TD+1); ++rr)
    {
        const VtxCoordType<GD>& XC_rr = VX->Get_Point(VI[rr]);
        for (SmallIndType cc=rr+1; cc < (TD+1); ++cc)
        {
            // compute the length of the current edge
            const VtxCoordType<GD>& XC_cc = VX->Get_Point(VI[cc]);
            const VtxCoordType<GD> DIFF = XC_rr.Subtract(XC_cc);
            const RealType LEN = DIFF.Length();
            if (LEN > MAX_Edge_Length) MAX_Edge_Length = LEN; // update
        }
    }
    return MAX_Edge_Length;
}

/***************************************************************************************/
/* compute the bounding (cartesian) box of a simplex */
template <SmallIndType TD, SmallIndType GD>
inline void Simplex_Bounding_Box(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                 PointType* BB_min, PointType* BB_max)
{
    const VtxIndType Num_Vtx = TD+1;
    VX.Bounding_Box(Num_Vtx, VI, BB_min, BB_max);
}

/***************************************************************************************/
/* compute the volume of a simplex */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Volume(const BasePtCoord<GD>* const& VX, const VtxIndType* VI)
{
    // get the affine map back to the reference simplex
    Eigen::Matrix<PointType, GD, TD> A;
    Eigen::Matrix<PointType, GD, 1>  b;
    Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
    
    RealType det_A = 0.0;
    if (TD==0)
    {
        // nothing to do!
    }
    else if (GD==TD)
    {
        det_A = A.determinant();
    }
    else
    {
        // must use the metric
        Eigen::Matrix<PointType, TD, TD> G = A.transpose() * A;
        det_A = std::sqrt( G.determinant() );
    }
    
    return det_A / (RealType) my_factorial(TD);
}

/***************************************************************************************/
/* compute the perimeter of a simplex */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Perimeter(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                  RealType Facet_Vol[TD+1])
{
    RealType Perimeter = 0.0;
    
    if (TD==0)
    {
        Facet_Vol[0] = 0.0;
    }
    else if (TD==1)
    {
        // counting measure
        Facet_Vol[0] = 1.0;
        Facet_Vol[1] = 1.0;
    }
    else
    {
        // loop through each facet of the simplex
        for (SmallIndType ff=0; ff < (TD+1); ++ff)
        {
            // get the vertex indices of the current facet
            VtxIndType Facet_VI[TD];
            SmallIndType Shift = 0;
            for (SmallIndType vv=0; vv < TD; ++vv)
            {
                if (vv==ff) Shift = 1;
                Facet_VI[vv] = VI[vv + Shift];
            }
            // compute it's volume (or "surface area")
            Facet_Vol[ff] = Simplex_Volume<TD-1,GD>(VX, Facet_VI);
            //std::cout << "Facet_Vol[ff]: " << Facet_Vol[ff] << std::endl;
        }
    }
    
    Perimeter = std::accumulate(Facet_Vol, Facet_Vol+(TD+1), 0.0);
    //std::cout << "Perimeter: " << Perimeter << std::endl;
    return Perimeter;
}

/***************************************************************************************/
/* compute the barycenter of a simplex */
template <SmallIndType TD, SmallIndType GD>
inline void Simplex_Barycenter(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                               PointType* PC)
{
    // all barycentric coordinates are the same
    const PointType BC_Coord = 1.0 / (PointType) (TD+1);

    // init to first contribution
    const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
    for (SmallIndType qq=0; qq < GD; ++qq)
        PC[qq] = BC_Coord * XC_0[gg];

    for (SmallIndType tt=1; tt < (TD+1); ++tt)
    {
        // get one of the point coordinates
        const PointType* XC_tt = VX->Get_Point_coord(VI[tt]);
        for (SmallIndType gg=0; gg < GD; ++gg)
        {
            // compute contribution
            PC[gg] += BC_Coord * XC_tt[gg];
        }
    }
}

/***************************************************************************************/
/* compute the circumcenter and circumradius of a simplex.
   see:  https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm#Coxeter
   Inputs: vertices of the simplex
   Outputs: barycentric coordinates of circumcenter, and circumradius */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Circumcenter(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                     PointType* CB)
{
    // allocate (Cayley-Menger) matrix
    Eigen::Matrix<RealType, TD+2, TD+2> MM;
    MM.setZero();

    // set the first row (first entry is zero)
    for (SmallIndType kk=1; kk < (TD+2); ++kk)
    {
        MM(0,kk) = 1.0;
    }

    // now set the rest, based on length of simplex edges (squared)
    for (SmallIndType rr=0; rr < (TD+1); ++rr)
    {
        // get the (rr)-th vertex
        const PointType* P_rr = VX->Get_Point_coord(VI[rr]);
        for (SmallIndType cc=rr+1; cc < (TD+1); ++cc)
        {
            // get the (cc)-th vertex
            const PointType* P_cc = VX->Get_Point_coord(VI[cc]);
            // take the difference
            Eigen::Matrix<RealType, GD, 1> DD;
            for (SmallIndType gg=0; gg < GD; ++gg)
            {
                DD(gg) = P_rr[gg] - P_cc[gg];
            }
            MM(rr+1,cc+1) = DD.squaredNorm(); // length squared
        }
    }

    // now set the lower triangular part (it is a symmetric matrix)
    for (SmallIndType rr=0; rr < (TD+2); ++rr)
    {
        for (SmallIndType cc=rr+1; cc < (TD+2); ++cc)
        {
            MM(cc,rr) = MM(rr,cc);
        }
    }
    //std::cout << "MM = " << std::endl;
    //std::cout << MM << std::endl;

    // compute the inverse
    Eigen::Matrix<RealType, TD+2, TD+2> MM_inv = MM.inverse();
    // error check
    if (MM_inv(0,0) > 0.0)
    {
        std::cerr << "Fatal error in 'Simplex_Circumcenter'!" << std::endl;
        std::cerr << "     The circumradius was invalid!" << std::endl;
        std::exit(1);
    }
    //std::cout << "MM_inv = " << std::endl;
    //std::cout << MM_inv << std::endl;

    // output barycentric coordinates
    for (SmallIndType tt=0; tt < (TD+1); ++tt)
    {
        CB[tt] = MM_inv(0,tt+1);
    }
    // output circumradius
    const RealType CR = std::sqrt(-MM_inv(0,0) / 2.0);
    return CR;
}

/***************************************************************************************/
/* compute the incenter and inradius of (simplex) cells in the mesh.
   see:  "Coincidences of simplex centers and related facial structures" by Edmonds, et al.
   Inputs: vertices of the simplex
   Outputs: barycentric coordinates of incenter, and inradius */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Incenter(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                 PointType* CB)
{
    RealType SA[TD+1]; // holds facet "surface areas"
    const RealType Total_Area = Simplex_Perimeter<TD,GD>(VX, VI, SA);

    if (TD==0)
    {
        // the incenter must coincide with the sole point
        CB[0] = 1.0;
    }
    else
    {
        // output barycentric coordinates
        for (SmallIndType tt=0; tt < (TD+1); ++tt)
        {
            // convert to barycentric
            SA[tt] = SA[tt] / Total_Area;
            CB[tt] = SA[tt];
        }
    }

    RealType CR = 0.0;
    if (TD==0)
    {
        // a point has zero radius
        CR = 0.0;
    }
    else if (TD==1)
    {
        // radius is half the length of the cell
        CR = (1.0/2.0) * Simplex_Volume<TD,GD>(VX, VI);
    }
    else
    {
        /* convert barycentric incenter to cartesian incenter */
        PointType IC_cart[GD];
        // init to 0th contribution
        const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
        for (SmallIndType qq=0; qq < GD; ++qq)
            IC_cart[qq] = SA[0] * XC_0[qq];

        for (SmallIndType tt=1; tt < (TD+1); ++tt)
        {
            // get one of the point coordinates
            const PointType* XC_tt = VX->Get_Point_coord(VI[tt]);
            for (SmallIndType gg=0; gg < GD; ++gg)
            {
                // compute contribution
                IC_cart[gg] += SA[tt] * XC_tt[gg];
            }
        }

        // compute distance from incenter to the 0th facet of the simplex (any facet would do)
        // this is the inradius:
        CR = Hyperplane_Closest_Point<TD-1,GD>(VX, VI, IC_cart);
    }
    // output inradius
    return CR;
}

/***************************************************************************************/
/* compute the "shape regularity" of the simplex,
           i.e. the ratio of the circumradius to the inradius.
   Inputs: vertices of the simplex
   Outputs: a number representing the "shape regularity" ratio. */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Shape_Regularity(const BasePtCoord<GD>* const& VX, const VtxIndType* VI)
{
    // hold the coordinates
    PointType CB[TD+1];
    
    const RealType Circumradius = Simplex_Circumcenter<TD,GD>(VX, VI, CB);
    const RealType Inradius     = Simplex_Incenter<TD,GD>(VX, VI, CB);
    
    // output inradius
    return (Circumradius / Inradius);
}

/***************************************************************************************/
/* get orthogonal column vectors that span \R^{GD}, where the first TD vectors span
   the tangent space of the simplex, and the last (GD - TD) vectors span the space
   *orthogonal* to the tangent space (i.e. the normal space).
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex.
   output: Eigen matrix, whose columns are the (unit) basis vectors. */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Orthogonal_Frame(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                         Eigen::Matrix<PointType, GD, GD>& Ortho)
{
    Eigen::Matrix<PointType, GD, TD> A;
    Eigen::Matrix<PointType, GD, 1>  b;
    Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
    
    Eigen::ColPivHouseholderQR<Eigen::Matrix<PointType, GD, TD>> QR(A);
    const RealType A_rank = QR.rank();
    
    // get the Q matrix
    Ortho = QR.householderQ();
    // make it have positive determinant
    const RealType O_det = Ortho.determinant();
    if (O_det < 0.0)
    {
        Ortho.col(0) = -Ortho.col(0);
    }
    
    // adjust for orientation when GD==TD+1, i.e. keep the orientation consistent with A.
    if (GD==TD+1)
    {
        Eigen::Matrix<PointType, GD, GD> Test_Mat;
        // append the normal vector to the Jacobian matrix
        Test_Mat << A, Ortho.col(GD-1);
        const RealType Test_det = Test_Mat.determinant();
        if (Test_det < 0)
        {
            // flip the orientation of the normal vector...
            Ortho.col(GD-1) = -Ortho.col(GD-1);
            // ... and the tangent space
            Ortho.col(0)    = -Ortho.col(0);
        }
    }
    // std::cout << "Ortho" << "\n\n";
    // std::cout << Ortho << "\n\n";
    // std::cout << "rank is " << A_rank << "\n\n";
    return A_rank;
}

/***************************************************************************************/
/* get orthogonal column vectors that span the tangent space of the simplex.
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex.
   output: Eigen matrix, whose columns are the (unit) basis vectors. */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Tangent_Space(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                      Eigen::Matrix<PointType, GD, TD>& TS)
{
    RealType TS_rank = 0.0; // init
    
    if (TD==0)
    {
        // do nothing; TS is empty
    }
    else if (TD==1)
    {
        // special case: 1-D space curve
        Eigen::Matrix<PointType, GD, TD> A;
        Eigen::Matrix<PointType, GD, 1>  b;
        Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
        
        // normalize the one column vector of A to get the tangent vector
        TS = A * (1.0 / A.norm());
        TS_rank = 1.0;
    }
    else
    {
        Eigen::Matrix<PointType, GD, GD> Ortho;
        TS_rank = Simplex_Orthogonal_Frame<TD,GD>(VX, VI, Ortho);
        
        TS.setIdentity();
        TS = Ortho * TS;
        std::cout << TS << "\n\n";
        std::cout << "rank is " << TS_rank << "\n\n";
    }
    return TS_rank;
}

/***************************************************************************************/
/* get orthogonal column vectors that span the space *orthogonal* to the tangent space
   of the simplex.
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex.
   output: Eigen matrix, whose columns are the (unit) basis vectors. */
template <SmallIndType TD, SmallIndType GD>
inline RealType Simplex_Normal_Space(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                     Eigen::Matrix<PointType, GD, (GD-TD)>& NS)
{
    RealType NS_rank = 0.0; // init
    
    if (TD==0)
    {
        // its all of \R^{GD}
        NS_rank = (RealType) GD;
        NS.setIdentity();
    }
    else if ( (TD==1) && (GD==2) )
    {
        // special case: curve in \R^2
        Eigen::Matrix<PointType, GD, TD> A;
        Eigen::Matrix<PointType, GD, 1>  b;
        Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
        // rotate the column vector from A
        Eigen::Vector2d V(-A(1),A(0));
        
        // normalize to get the normal vector
        NS = V * (1.0 / V.norm());
        NS_rank = 1.0;
    }
    else if ( (TD==2) && (GD==3) )
    {
        // special case: surface in \R^3
        Eigen::Matrix<PointType, GD, TD> A;
        Eigen::Matrix<PointType, GD, 1> b;
        Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
        // take cross product of the two column vectors of A
        Eigen::Matrix<PointType, GD, 1> V; // = A.col(0).cross(A.col(1));
        // must do this manually because of templating issue
        V(0) =   A(1,0) * A(2,1) - A(2,0) * A(1,1);
        V(1) = -(A(0,0) * A(2,1) - A(2,0) * A(0,1));
        V(2) =   A(0,0) * A(1,1) - A(1,0) * A(0,1);
        
        // normalize to get the normal vector
        NS = V * (1.0 / V.norm());
        NS_rank = 1.0;
    }
    else
    {
        // general case
        Eigen::Matrix<PointType, GD, GD> Ortho;
        NS_rank = Simplex_Orthogonal_Frame<TD,GD>(VX, VI, Ortho);
        if (TD==GD) NS_rank = 0.0;
        
        // get the last columns
        for (SmallIndType jj=0; jj < (GD-TD); ++jj)
        {
            for (SmallIndType ii=0; ii < GD; ++ii)
            {
                NS(ii,jj) = Ortho(ii,TD+jj);
            }
            // if (GD > TD) NS.col(jj) = Ortho.col(TD+jj); // this doesn't work
        }
        std::cout << NS << "\n\n";
        std::cout << "rank is " << NS_rank << "\n\n";
    }
    return NS_rank;
}

// /***************************************************************************************/
// /* get the angles that the facets of the simplex make with respect to each other.
   // inputs: global list of vertex coordinates, ordered list of vertices of the simplex.
   // output: an array of angles of length ??????.
   // NOTE: this only works when TD==GD; the rest is to be implemented.  */
// template <SmallIndType TD, SmallIndType GD>
// void Simplex_Angles(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                    // RealType* Angles)
// {
    // RealType NS_rank = 0.0; // init
    
    // if (TD==0)
    // {
        // // its all of \R^{GD}
        // NS_rank = (RealType) GD;
        // NS.setIdentity();
    // }
    // else if ( (TD==1) && (GD==2) )
    // {
        // // special case: curve in \R^2
        // Eigen::Matrix<PointType, GD, TD> A;
        // Eigen::Matrix<PointType, GD, 1>  b;
        // Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
        // // rotate the column vector from A
        // Eigen::Vector2d V(-A(1),A(0));
        
        // // normalize to get the normal vector
        // NS = V * (1.0 / V.norm());
        // NS_rank = 1.0;
    // }
    // else if ( (TD==2) && (GD==3) )
    // {
        // // special case: surface in \R^3
        // Eigen::Matrix<PointType, GD, TD> A;
        // Eigen::Matrix<PointType, GD, 1> b;
        // Simplex_Affine_Map<TD,GD>(VX, VI, A, b);
        // // take cross product of the two column vectors of A
        // Eigen::Matrix<PointType, GD, 1> V; // = A.col(0).cross(A.col(1));
        // // must do this manually because of templating issue
        // V(0) =   A(1,0) * A(2,1) - A(2,0) * A(1,1);
        // V(1) = -(A(0,0) * A(2,1) - A(2,0) * A(0,1));
        // V(2) =   A(0,0) * A(1,1) - A(1,0) * A(0,1);
        
        // // normalize to get the normal vector
        // NS = V * (1.0 / V.norm());
        // NS_rank = 1.0;
    // }
    // else
    // {
        // // general case
        // Eigen::Matrix<PointType, GD, GD> Ortho;
        // NS_rank = Simplex_Orthogonal_Frame<TD,GD>(VX, VI, Ortho);
        // if (TD==GD) NS_rank = 0.0;
        
        // // get the last columns
        // for (SmallIndType jj=0; jj < (GD-TD); ++jj)
        // {
            // for (SmallIndType ii=0; ii < GD; ++ii)
            // {
                // NS(ii,jj) = Ortho(ii,TD+jj);
            // }
            // // if (GD > TD) NS.col(jj) = Ortho.col(TD+jj); // this doesn't work
        // }
        // std::cout << NS << "\n\n";
        // std::cout << "rank is " << NS_rank << "\n\n";
    // }
    // return NS_rank;
// }

/***************************************************************************************/
/* find closest point, X_star, on a hyperplane to a given point, Y.
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex
           that defines the hyperplane, the given point (cartesian) coordinates (as an array).
   output: Eigen::Vector of the projection of (X_0 - Y) onto normal space of simplex;
           X_0 is the 0th vertex of the simplex.
   Note: the closest point may not lie inside the simplex.  */
template <SmallIndType TD, SmallIndType GD>
inline void Hyperplane_Closest_Point(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                     const PointType* PY, Eigen::Matrix<PointType, GD, 1>& PROJ)
{
    PROJ.setZero(); // init
    
    if (GD==TD)
    {
        // easy case: it's the same as PY
        // so PROJ is the zero vector
    }
    else
    {
        Eigen::Matrix<PointType, GD, (GD-TD)> NS;
        Simplex_Normal_Space<TD,GD>(VX, VI, NS);
        
        // get one of the points of the simplex
        const PointType* XC_0 = VX->Get_Point_coord(VI[0]);
        
        // compute the difference vector
        Eigen::Matrix<PointType, GD, 1> DIFF;
        for (SmallIndType kk=0; kk < GD; ++kk)
        {
            DIFF(kk) = XC_0[kk] - PY[kk];
        }
        
        // project onto the normal space
        for (SmallIndType qq=0; qq < (GD-TD); ++qq)
        {
            PROJ += NS.col(qq).dot(DIFF) * NS.col(qq);
        }
    }
}
/* find closest point, X_star, on a hyperplane to a given point, Y.
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex
           that defines the hyperplane, the given point (cartesian) coordinates (as an array).
   output: closest point on hyperplane in cartesian coordinates (as an array);
           the distance between X_star and Y is also returned.
   Note: the closest point may not lie inside the simplex.  */
template <SmallIndType TD, SmallIndType GD>
inline RealType Hyperplane_Closest_Point(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                         const PointType* PY, PointType* PX_star)
{
    Eigen::Matrix<PointType, GD, 1> PROJ;
    Hyperplane_Closest_Point<TD,GD>(VX, VI, PY, PROJ); // PROJ = P(X_0 - Y)
    RealType DIST = PROJ.norm();
    
    // get the closest point
    for (SmallIndType kk=0; kk < GD; ++kk)
    {
        PX_star[kk] = PROJ(kk) + PY[kk];
    }
    
    return DIST;
}
/* find minimum distance between a hyperplane and a given point, Y.
   inputs: global list of vertex coordinates, ordered list of vertices of the simplex
           that defines the hyperplane, the given point (cartesian) coordinates (as an array).
   output: minimum distance between hyperplane and Y.
   Note: the closest point may not lie inside the simplex.  */
template <SmallIndType TD, SmallIndType GD>
inline RealType Hyperplane_Closest_Point(const BasePtCoord<GD>* const& VX, const VtxIndType* VI,
                                         const PointType* PY)
{
    Eigen::Matrix<PointType, GD, 1> PROJ;
    Hyperplane_Closest_Point<TD,GD>(VX, VI, PY, PROJ); // PROJ = P(X_0 - Y)
    const RealType DIST = PROJ.norm();
    return DIST;
}





/***/
