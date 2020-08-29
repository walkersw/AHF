/*
============================================================================================
   Some simple mathematical sub-routines.

   Copyright (c) 07-15-2020,  Shawn W. Walker
============================================================================================
*/

// std C++ classes
//#include <algorithm>
//#include <vector>
//#include <map>
//#include <numeric>

#define _BASICMATH_H

#ifndef _PRELIM_H
#include "Prelim.h"  // simple typedefs, etc...
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

/***************************************************************************************/
/* saturation for dot product */
inline RealType angle_dot_prod_normalize(RealType& DP)
{
    if (DP < -1.0)
        return -1.0;
    else if (DP > 1.0)
        return 1.0;
    else
        return DP;
}

// local to global mapping...  put it into a different file.


/***/
