/*
============================================================================================
   This is a scrap unit test for the AHF data structure implementation.

   This example...

   Copyright (c) 11-20-2015,  Shawn W. Walker
============================================================================================
*/

#include "src_code/TriMesh.cc"

using namespace std;

// test it!
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, -1 is failure

    // create the object
    TriMesh  TM;

    BaseMesh<1> TEST1;
    BaseMesh<3> TEST3;
    BaseMesh<4> TEST4;

    // start over and make a non-manifold mesh
    TM.Clear();

    // make a long list of cells to test speed
    const CellIndType MAX = 20000;
    TM.Reserve_Cell(MAX);
    for (CellIndType ii = 1; ii <= MAX; ii++)
    {
        TM.Append_Cell(3*(ii-1) + 1, 3*(ii-1) + 2, 3*(ii-1) + 3);
    }
    TM.Finalize_v2hfs();

    if (OUTPUT_CODE==0)
        cout << "Unit test is successful!" << endl;
    else
        cout << "Unit test failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
