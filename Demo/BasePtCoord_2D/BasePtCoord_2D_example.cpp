/*
============================================================================================
   This is a demo for how to use the BasePtCoord class.  It is just a glorified
   C++ STL vector.

   This example stores some points (in 2-D) and does some basic processing.  The usage
   for other dimensions is obvious.

   Note: see "BasePtCoord.cc" comments for more info.

   Copyright (c) 02-24-2020,  Shawn W. Walker
============================================================================================
*/

#include "../../src_code/TypedefMeshes.h"

using namespace std;

// demo
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, > 0 is failure

    // create the object: a set of point coordinates in 2-D
    BasePtCoord<2>  BP;

    // define the points (5 points)
    BP.Init_Points(5); // initialize all 5 points to the origin
    BP.Set_Coord(0,-1.3,-2.4); // set 0th point to (-1.3,-2.4)
    BP.Set_Coord(1, 0.9,-0.8); // set 1st point to ( 0.9,-0.8)
    BP.Set_Coord(2, 0.5, 1.1); // etc...
    BP.Set_Coord(3,-1.6, 1.4);
    BP.Set_Coord(4, 0.1,-0.04);
	cout << endl;

	// now close it to further modification
	BP.Close();

    // let's not risk changing the points inadvertently
    const BasePtCoord<2>& BP_c = BP;

    // basic properties
    const SmallIndType GD = BP_c.Geo_Dim();
    const VtxIndType NP = BP_c.Num_Points();
    if (GD!=2)
    {
        cout << "The dimension of the points does not equal 2.  Error!" << endl;
        OUTPUT_CODE = 1;
    }
    if (NP!=5)
    {
        cout << "The number of points does not equal 5.  Error!" << endl;
        OUTPUT_CODE = 2;
    }
    
    // access the 2nd point's coordinates
    // NOTE: "PointType" is just a double, but this library is really anal about types!
    const VtxIndType my_pt = 2;
    const PointType* PC = BP_c.Get_Point_coord(my_pt);
	cout << "The coordinates of Point #" << my_pt << ": ";
    cout << "[" << PC[0] << ", " << PC[1] << "]" << endl;
    // define an object to hold this point
    VtxCoordType<2> PC_VTX;
    PC_VTX.Set(PC);
    cout << endl;
    
    // access the 2nd point another way
    const VtxCoordType<2> VV = BP_c.Get_Point(my_pt);
	cout << "Point #" << my_pt << " is described by: " << endl;
    VV.Print();
    cout << endl;
    
    // can also output the coordinates "manually"
	cout << "Point #" << my_pt << " has the following coordinates (manual):" << endl;
    cout << "[" << VV.coord[0] << ", " << VV.coord[1] << "]" << endl;
	cout << endl;
    
    // do a comparison
    if (VV.Equal(PC_VTX))
    {
        cout << "Both versions of the 2nd point are equal." << endl;
    }
    else
    {
        cout << "The two versions of the 2nd point are NOT equal!" << endl;
        OUTPUT_CODE = 3;
    }
    cout << endl;

    // print them all out
    BP_c.Display_Vtx_Coord();
    cout << endl;

    // now clear it all
    BP.Clear(); // use the original object here
    if (BP_c.Num_Points()==0)
    {
        cout << "Points are now empty." << endl;
    }
    else
    {
        cout << "The points are supposed to be empty, but they are NOT!" << endl;
        OUTPUT_CODE = 4;
    }
    cout << endl;

    if (OUTPUT_CODE==0)
        cout << "Demo completed successfully!" << endl;
    else
        cout << "Demo failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
