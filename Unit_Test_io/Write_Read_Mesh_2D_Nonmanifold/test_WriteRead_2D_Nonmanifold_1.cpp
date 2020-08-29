/*
============================================================================================
   This is a unit test for the Mesh_IO_VTK class.

   This example writes a 2-D non-manifold surface mesh to an ASCII/BINARY VTK file,
   and reads it back.

   Copyright (c) 08-28-2020,  Shawn W. Walker
============================================================================================
*/

#include "../src_code/TypedefMeshes.h"
#include "../io_code/io_vtk.cc"

using namespace std;

// unit test
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, 1 is failure

    // create the object: non-manifold 2-D surface mesh example
    BasePtCoord<3> VX;
    SurfaceMesh(TM,VX);
    //Mesh<2,3>  TM(&VX); // alternative

    // non-manifold 2-D mesh example
    TM.Reserve_Cells(4);
    TM.Append_Cell(0,1,2);
    TM.Append_Cell(1,3,2);
    TM.Append_Cell(1,4,2);
    TM.Append_Cell(1,2,5);

    // now add the vertex point coordinates
    VX.Init_Points(6);
    VX.Set_Coord(0, -1.0, 0.5, 0.0);
    VX.Set_Coord(1,  0.0, 0.0, 0.0);
    VX.Set_Coord(2,  0.0, 1.0, 0.0);
    VX.Set_Coord(3,  1.0, 0.5, 0.0);
    VX.Set_Coord(4,  0.0, 0.5, 1.0);
    VX.Set_Coord(5,  0.0, 0.5,-1.0);

    // create mesh writer
    Mesh_IO_VTK  MeshIO;

    // write the mesh data
    MeshIO.Write_ASCII_VTK("test_non_manifold_mesh_ASCII.vtk", "simple non-manifold mesh", TM, VX);
    MeshIO.Write_Binary_VTK("test_non_manifold_mesh_BIN.vtk", "simple non-manifold mesh", TM, VX);

    // read the mesh data ASCII (read into temporary structures)
    BasePtCoord<3> VX_Read_ASCII;
    SurfaceMesh(TM_Read_ASCII,VX_Read_ASCII);
    MeshIO.Read_ASCII_VTK("test_non_manifold_mesh_ASCII.vtk", TM_Read_ASCII);

    // check that cells are the same!
    SmallIndType TD = TM_Read_ASCII.Top_Dim();
    if (TD!=2)
    {
        cout << "(From ASCII File) Topological Dimension: " << TD << endl;
        cout << "                  It should be 2!" << endl;
        OUTPUT_CODE = 1;
    }
    CellIndType  NC = TM_Read_ASCII.Num_Cells();
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        //cout << "Check Cell #" << ci << endl;
        const CellType& CL = TM.Get_Cell(ci);
        const CellType& CL_ASCII = TM_Read_ASCII.Get_Cell(ci);
        if (!(CL.Equal(CL_ASCII)))
        {
            cout << "ASCII Check:" << endl;
            cout << "Cell #" << ci << " does not match!" << endl;
            CL.Print();
            CL_ASCII.Print();
            OUTPUT_CODE = 2;
            break;
        }
    }
    // check that vertex coordinates are the same!
    VtxIndType NP = VX_Read_ASCII.Num_Points();
    for (VtxIndType vi = 0; vi < NP; ++vi)
    {
        //cout << "Vtx #" << vi << endl;
        const CoordType& VC = VX.Get_Point(vi);
        const CoordType& VC_ASCII = VX_Read_ASCII.Get_Point(vi);
        if (!(VC.Equal(VC_ASCII)))
        {
            cout << "ASCII Check:" << endl;
            cout << "Vtx #" << vi << " coordinates do not match!" << endl;
            VC.Print();
            VC_ASCII.Print();
            OUTPUT_CODE = 3;
            break;
        }
    }

    // read the mesh data BINARY (read into temporary structures)
    BasePtCoord<3> VX_Read_BIN;
    SurfaceMesh(TM_Read_BIN,VX_Read_BIN);
    MeshIO.Read_BINARY_VTK("test_non_manifold_mesh_BIN.vtk", TM_Read_BIN);

    // check that cells are the same!
    TD = TM_Read_BIN.Top_Dim();
    if (TD!=2)
    {
        cout << "(From BINARY File) Topological Dimension: " << TD << endl;
        cout << "                   It should be 2!" << endl;
        OUTPUT_CODE = 4;
    }
    NC = TM_Read_BIN.Num_Cells();
    for (CellIndType ci = 0; ci < NC; ++ci)
    {
        //cout << "Check Cell #" << ci << endl;
        const CellType& CL = TM.Get_Cell(ci);
        const CellType& CL_BIN = TM_Read_BIN.Get_Cell(ci);
        if (!(CL.Equal(CL_BIN)))
        {
            cout << "BINARY Check:" << endl;
            cout << "Cell #" << ci << " does not match!" << endl;
            CL.Print();
            CL_BIN.Print();
            OUTPUT_CODE = 5;
            break;
        }
    }
    // check that vertex coordinates are the same!
    NP = VX_Read_BIN.Num_Points();
    for (VtxIndType vi = 0; vi < NP; ++vi)
    {
        //cout << "Vtx #" << vi << endl;
        const CoordType& VC = VX.Get_Point(vi);
        const CoordType& VC_BIN = VX_Read_BIN.Get_Point(vi);
        if (!(VC.Equal(VC_BIN)))
        {
            cout << "BINARY Check:" << endl;
            cout << "Vtx #" << vi << " coordinates do not match!" << endl;
            VC.Print();
            VC_BIN.Print();
            OUTPUT_CODE = 6;
            break;
        }
    }

    if (OUTPUT_CODE==0)
        cout << "Unit test is successful!" << endl;
    else
        cout << "Unit test failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
