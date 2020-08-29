/*
============================================================================================
   This is a unit test for the AHF data structure implementation.

   This example stores a 2-D manifold mesh (consisting of 4 triangles) that has
   0 non-manifold edges and 0 non-manifold vertices.  It also does some basic processing.
   See 'Manifold_Triangle_Mesh_1.jpg' for a picture of the mesh embedded in 2-D.

   Copyright (c) 05-21-2020,  Shawn W. Walker
============================================================================================
*/

#include <TypedefMeshes.h>

using namespace std;

// unit test
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, > 0 is failure

    // create the object: manifold 2-D mesh example
    BasePtCoord<2> VX;
    TriMesh(TM,VX);
    //Mesh<2,2>  TM(&VX); // alternative
    
    // define the cell connectivity
    TM.Reserve_Cells(4);
    TM.Append_Cell(0,1,4);
    TM.Append_Cell(1,2,4);
    TM.Append_Cell(2,3,4);
    TM.Append_Cell(3,0,4);

    // now add the vertex point coordinates
    VX.Init_Points(5);
    VX.Set_Coord(0, 0.0,0.0);
    VX.Set_Coord(1, 1.0,0.0);
    VX.Set_Coord(2, 1.0,1.0);
    VX.Set_Coord(3, 0.0,1.0);
    VX.Set_Coord(4, 0.5,0.5);

    // // re-index
    // const VtxIndType old_indices[5] = {0, 1, 2, 3, 4};
    // const VtxIndType new_indices[5] = {3, 7, 8, 9, 12};
    // TM.Reindex_Vertices(5, new_indices);
    // VX.Reindex_Vertices(5, old_indices, new_indices);

    // now display coordinates
    cout << endl;
    VX.Display_Vtx_Coord();

    // we now stop adding cells and coordinates
    TM.Finalize_v2hfs_DEBUG();
    
    // now display the half-facets (attached to vertices) of the intermediate structure 'v2hfs'
    cout << endl;
    TM.Display_v2hfs();

    // check 'v2hfs' against reference data
    VtxHalfFacetType v2hfs_REF[12];
    v2hfs_REF[ 0].Set(1,0,2);
    v2hfs_REF[ 1].Set(2,1,2);
    v2hfs_REF[ 2].Set(3,2,2);
    v2hfs_REF[ 3].Set(3,3,2);
    v2hfs_REF[ 4].Set(4,0,0);
    v2hfs_REF[ 5].Set(4,0,1);
    v2hfs_REF[ 6].Set(4,1,0);
    v2hfs_REF[ 7].Set(4,1,1);
    v2hfs_REF[ 8].Set(4,2,0);
    v2hfs_REF[ 9].Set(4,2,1);
    v2hfs_REF[10].Set(4,3,0);
    v2hfs_REF[11].Set(4,3,1);
    // error check
    const Vtx2HalfFacet_Mapping& c_V2HF_Map = TM.Get_v2hfs();
    const std::vector<VtxHalfFacetType>& c_v2hfs = c_V2HF_Map.Get_VtxMap();
    for (VtxIndType jj = 0; jj < c_v2hfs.size(); ++jj)
        if (!c_v2hfs[jj].Equal(v2hfs_REF[jj]))
        {
            cout << "Intermediate data 'v2hfs' failed!" << endl;
            OUTPUT_CODE = 1;
            break;
        }

    // fill out the sibling half-facet data in 'Cell'
    TM.Build_Sibling_HalfFacets_DEBUG();

    // display that cell connectivity data
    cout << endl;
    TM.Display_Cell();
    cout << endl;

    // check 'Cell' against reference data
    CellSimplexType<2> Cell_REF[4];
    HalfFacetType hf;
    // Cell #1
    Cell_REF[ 0].Set(0,0,hf.Set(1,1));
    Cell_REF[ 0].Set(1,1,hf.Set(3,0));
    Cell_REF[ 0].Set(2,4,hf.Set());
    // Cell #2
    Cell_REF[ 1].Set(0,1,hf.Set(2,1));
    Cell_REF[ 1].Set(1,2,hf.Set(0,0));
    Cell_REF[ 1].Set(2,4,hf.Set());
    // Cell #3
    Cell_REF[ 2].Set(0,2,hf.Set(3,1));
    Cell_REF[ 2].Set(1,3,hf.Set(1,0));
    Cell_REF[ 2].Set(2,4,hf.Set());
    // Cell #4
    Cell_REF[ 3].Set(0,3,hf.Set(0,1));
    Cell_REF[ 3].Set(1,0,hf.Set(2,0));
    Cell_REF[ 3].Set(2,4,hf.Set());

    // error check
    for (CellIndType cc = 0; cc < TM.Num_Cells(); ++cc)
    {
        const VtxIndType* C_vtx = NULL;
        const HalfFacetType* C_hf = NULL;
        TM.Get_Cell(cc, C_vtx, C_hf);
        if (!Cell_REF[cc].Equal(C_vtx, C_hf))
        {
            cout << "Cell connectivity (and sibling half-facets) is incorrect!" << endl;
            OUTPUT_CODE = 2;
            break;
        }
    }

    // generate the Vtx2HalfFacets data struct
    TM.Build_Vtx2HalfFacets_DEBUG();
    TM.Close(); // we can close it now, i.e. no more modifications

    // now display half-facets attached to vertices for final data structure
    TM.Display_Vtx2HalfFacets();
    cout << endl;

    // check 'Vtx2HalfFacets' against reference data
    VtxHalfFacetType Vtx2HalfFacets_REF[5];
    Vtx2HalfFacets_REF[ 0].Set( 0,3,2);
    Vtx2HalfFacets_REF[ 1].Set( 1,1,2);
    Vtx2HalfFacets_REF[ 2].Set( 2,2,2);
    Vtx2HalfFacets_REF[ 3].Set( 3,3,2);
    Vtx2HalfFacets_REF[ 4].Set( 4,0,0);
    // error check
    const std::vector<VtxHalfFacetType>& c_Vtx2HF = TM.Get_Vtx2HalfFacets().Get_VtxMap();
    for (VtxIndType jj = 0; jj < c_Vtx2HF.size(); ++jj)
    {
        if (!c_Vtx2HF[jj].Equal(Vtx2HalfFacets_REF[jj]))
        {
            cout << "'Vtx2HalfFacets' data is incorrect!" << endl;
            OUTPUT_CODE = 3;
            break;
        }
    }

    // display cells attached to a vertex
    VtxIndType V_IN = 2;
    TM.Display_Cells_Attached_To_Vertex(V_IN);
    cout << endl;

    // check attached cells against reference data
    std::vector<CellIndType> cell_ind_1;
    TM.Get_Cells_Attached_To_Vertex(V_IN, 2, cell_ind_1);
    const CellIndType cell_ind_1_REF[2] = {2, 1};
    for (unsigned int jj = 0; jj < cell_ind_1.size(); ++jj)
        if (cell_ind_1[jj]!=cell_ind_1_REF[jj])
        {
            cout << "Cell attachment data is incorrect!" << endl;
            OUTPUT_CODE = 4;
            break;
        }

    // display if two cells are facet connected
    const  VtxIndType V1 = 1;
    const CellIndType C1 = 0;
    const CellIndType C2 = 1;
    TM.Display_Two_Cells_Are_Facet_Connected(V1, C1, C2);
    cout << endl;

    // check facet connected cells against reference data
    const bool CHK_Facet_Connected = TM.Two_Cells_Are_Facet_Connected(V1, C1, C2);
    const bool CHK_Facet_Connected_REF = true;
    if (CHK_Facet_Connected!=CHK_Facet_Connected_REF)
    {
        cout << "Facet connected cell data is incorrect!" << endl;
        OUTPUT_CODE = 5;
    }

    // display half-facets attached to given half-facet
    HalfFacetType TEST_attached;
    TEST_attached.Set(3,1);
    TM.Display_HalfFacets_Attached_To_HalfFacet(TEST_attached);
    cout << endl;

    // check attached half-facets against reference data
    HalfFacetType attached_REF[2];
    attached_REF[0].Set(3,1);
    attached_REF[1].Set(2,0);
    std::vector<HalfFacetType> attached;
    TM.Get_HalfFacets_Attached_To_HalfFacet(TEST_attached, attached);
    for (unsigned int jj = 0; jj < attached.size(); ++jj)
        if (!attached[jj].Equal(attached_REF[jj]))
        {
            cout << "HalfFacet attachment data is incorrect!" << endl;
            OUTPUT_CODE = 6;
            break;
        }

    // check that half-facet has no neighbor
    TEST_attached.Set(2,2);
    TM.Get_HalfFacets_Attached_To_HalfFacet(TEST_attached, attached);
    if (attached.size()!=1)
    {
        cout << "HalfFacet attachment data is incorrect!" << endl;
        OUTPUT_CODE = 7;
    }
    else if (!TEST_attached.Equal(attached[0]))
    {
        cout << "HalfFacet attachment should only include the initial given half-facet..." << endl;
        cout << "          but it does not!" << endl;
        OUTPUT_CODE = 8;
    }

    // display non-manifold facets (edges)
    TM.Display_Nonmanifold_HalfFacets();
    cout << endl;

    // check non-manifold facets (edges) against reference data
    std::vector<HalfFacetType> non_manifold_hf;
    TM.Get_Nonmanifold_HalfFacets(non_manifold_hf);
    if (non_manifold_hf.size()!=0)
    {
        cout << "Non-manifold facet data is incorrect!" << endl;
        OUTPUT_CODE = 9;
    }

    // display non-manifold vertices
    TM.Display_Nonmanifold_Vertices();
    cout << endl;

    // check non-manifold vertices against reference data
    std::vector<VtxIndType> non_manifold_vtx;
    TM.Get_Nonmanifold_Vertices(non_manifold_vtx);
    if (non_manifold_vtx.size()!=0)
    {
        cout << "Non-manifold vertex data is incorrect!" << endl;
        OUTPUT_CODE = 10;
    }

    TM.Display_Unique_Vertices();
    cout << endl;

    // check unique vertices against reference data
    std::vector<VtxIndType> uv;
    TM.Get_Unique_Vertices(uv);
    const VtxIndType uv_REF[5] = {0, 1, 2, 3, 4};
    for (unsigned int jj = 0; jj < uv.size(); ++jj)
        if (uv[jj]!=uv_REF[jj])
        {
            cout << "Unique vertex data is incorrect!" << endl;
            OUTPUT_CODE = 11;
            break;
        }

    // display all edges in the mesh
    std::vector<MeshEdgeType> edges;
    TM.Get_Edges(edges);
    cout << "A unique list of edges in the mesh:" << endl;
    std::vector<MeshEdgeType>::const_iterator it;
    for (it = edges.begin(); it!=edges.end(); ++it)
    {
        cout << "[" << (*it).vtx[0] << ", " << (*it).vtx[1] << "]" << endl;
    }
    cout << endl;

    // test if a pair of vertices is connected by an edge
    cout << "Vtx #0 and Vtx #4 are connected." << endl;
    if (!TM.Is_Connected(0, 4))
    {
        cout << "Incorrect!  They are NOT connected!" << endl;
        OUTPUT_CODE = 12;
    }
    // test if a pair of vertices is connected by an edge
    cout << "Vtx #0 and Vtx #2 are NOT connected." << endl;
    if (TM.Is_Connected(0, 2))
    {
        cout << "Incorrect!  They ARE connected!" << endl;
        OUTPUT_CODE = 13;
    }
    cout << endl;

    // display all cells attached to an edge
    std::vector<CellIndType> attached_cells;
    MeshEdgeType EE;
    EE.Set(1,4);
    TM.Get_Cells_Attached_To_Edge(EE, attached_cells);
    cout << "Here are the cell indices of cells attached to the edge [1, 4]:" << endl;
    std::vector<CellIndType>::const_iterator ic;
    for (ic = attached_cells.begin(); ic!=attached_cells.end(); ++ic)
        cout << "Cell #" << (*ic) << endl;
    cout << endl;

    // display the free boundary
    std::vector<HalfFacetType> free_bdy;
    TM.Get_FreeBoundary(free_bdy);
    cout << "Here are the half-facets that lie on the free boundary:" << endl;
    std::vector<HalfFacetType>::const_iterator hfi;
    for (hfi = free_bdy.begin(); hfi!=free_bdy.end(); ++hfi)
    {
        (*hfi).Print();
        cout << endl;
    }
    cout << endl;
    
    // test: Barycentric_To_Cartesian
    const CellIndType CI_A[2] = {0, 2};
    const PointType   PB_in[6] = {(1.0/3.0), (1.0/3.0), (1.0/3.0),   0.2, 0.6, 0.2}; // two sets of barycentric coordinates
    PointType  PC_out[4];
    TM.Barycentric_To_Cartesian(2, CI_A, PB_in, PC_out);
    cout << "These are the two points (from given barycenters):" << endl;
    for (CellIndType ii = 0; ii < 2; ++ii)
        cout << "Point (cartesian): (" << PC_out[2*ii + 0] << ", " << PC_out[2*ii + 1] << ")" << endl;
    const bool PC_CHK = (abs(PC_out[0]-0.5) > 1E-14) || (abs(PC_out[1]-(1.0/6.0)) > 1E-14)
                     || (abs(PC_out[2]-0.3) > 1E-14) || (abs(PC_out[3]-0.9) > 1E-14);
    if (PC_CHK)
    {
        cout << "Incorrect!  The two points should be: (0.5, 0.1666667) and (0.3, 0.9)!" << endl;
        OUTPUT_CODE = 14;
    }
    cout << endl;
    
    // test: Cartesian_To_Barycentric
    const CellIndType CI_B[2] = {1, 3};
    const PointType   PC_in[4] = {0.75, 0.5,   0.2, 0.6}; // two sets of cartesian coordinates
    PointType  PB_out[6];
    TM.Cartesian_To_Barycentric(2, CI_B, PC_in, PB_out);
    cout << "These are the two points (from given cartesian coordinates):" << endl;
    for (CellIndType ii = 0; ii < 2; ++ii)
        cout << "Point (barycentric): (" << PB_out[3*ii + 0] << ", " << PB_out[3*ii + 1] << ", " << PB_out[3*ii + 2] << ")" << endl;
    const bool PB_CHK = (abs(PB_out[0]-0.25) > 1E-14) || (abs(PB_out[1]-0.25) > 1E-14) || (abs(PB_out[2]-0.5) > 1E-14)
                     || (abs(PB_out[3]-0.4)  > 1E-14) || (abs(PB_out[4]-0.2)  > 1E-14) || (abs(PB_out[5]-0.4) > 1E-14);
    if (PB_CHK)
    {
        cout << "Incorrect!  The two points should be: (0.25, 0.25, 0.5) and (0.4, 0.2, 0.4)!" << endl;
        OUTPUT_CODE = 15;
    }
    cout << endl;
    
    // test: Circumcenter
    PointType  CB_out[6];
    RealType   CR_out[2];
    TM.Circumcenter(2, CI_A, CB_out, CR_out);
    cout << "Here are circumcenters and circumradii for cell #s: " << CI_A[0] << " and " << CI_A[1] << ":" << endl;
    for (CellIndType ii = 0; ii < 2; ++ii)
        cout << "circumcenter (barycentric): (" << CB_out[3*ii + 0] << ", " << CB_out[3*ii + 1] << ", " << CB_out[3*ii + 2] << ")" 
             << ", circumradius: " << CR_out[ii] << endl;
    const bool CB_CHK = (abs(CB_out[0]-0.5) > 1E-14) || (abs(CB_out[1]-0.5) > 1E-14) || (abs(CB_out[2]-0.0) > 1E-14)
                     || (abs(CB_out[3]-0.5) > 1E-14) || (abs(CB_out[4]-0.5) > 1E-14) || (abs(CB_out[5]-0.0) > 1E-14)
                     || (abs(CR_out[0]-0.5) > 1E-14) || (abs(CR_out[1]-0.5) > 1E-14);
    if (CB_CHK)
    {
        cout << "Incorrect!  The two circumcenters should be: (0.5, 0.5, 0.0) and (0.5, 0.5, 0.0)," << endl;
        cout << "            with a circumradius of 0.5 each!" << endl;
        OUTPUT_CODE = 16;
    }
    cout << endl;
    
    // test: incenter
    PointType  IB_out[6];
    RealType   IR_out[2];
    TM.Incenter(2, CI_A, IB_out, IR_out);
    cout << "Here are incenters and inradii for cell #s: " << CI_A[0] << " and " << CI_A[1] << ":" << endl;
    for (CellIndType ii = 0; ii < 2; ++ii)
        cout << "incenter (barycentric): (" << IB_out[3*ii + 0] << ", " << IB_out[3*ii + 1] << ", " << IB_out[3*ii + 2] << ")" 
             << ", inradius: " << IR_out[ii] << endl;
    const bool IB_CHK = (abs(IB_out[0]-0.292893218813452) > 1E-14) || (abs(IB_out[1]-0.292893218813452) > 1E-14)
                     || (abs(IB_out[2]-0.414213562373095) > 1E-14)
                     || (abs(IB_out[3]-0.292893218813452) > 1E-14) || (abs(IB_out[4]-0.292893218813452) > 1E-14)
                     || (abs(IB_out[5]-0.414213562373095) > 1E-14)
                     || (abs(IR_out[0]-0.207106781186548) > 1E-14) || (abs(IR_out[1]-0.207106781186548) > 1E-14);
    if (IB_CHK)
    {
        cout << "Incorrect!  The two incenters should be (both): ( (1/sqrt(2))/(1+sqrt(2)), (1/sqrt(2))/(1+sqrt(2)), 1/(1+sqrt(2)))," << endl;
        cout << "            with an inradius of (1/2) / (1+sqrt(2)) for both!" << endl;
        OUTPUT_CODE = 17;
    }
    cout << endl;
    
    // test: diameter
    const CellIndType CI_all[4] = {0, 1, 2, 3};
    RealType   Diam_out[4];
    TM.Diameter(4, CI_all, Diam_out);
    cout << "Here are the diameters for these cells: " << endl;
    for (CellIndType ii = 0; ii < 4; ++ii)
        cout << "Cell #" << CI_all[ii] << ",  Diameter: " << Diam_out[ii] << endl;
    const bool DIAM_CHK = (abs(Diam_out[0]-1.0) > 1E-14) || (abs(Diam_out[1]-1.0) > 1E-14)
                       || (abs(Diam_out[2]-1.0) > 1E-14) || (abs(Diam_out[3]-1.0) > 1E-14);
    if (DIAM_CHK)
    {
        cout << "Incorrect!  The diameters should be 1.0 for all." << endl;
        OUTPUT_CODE = 18;
    }
    cout << endl;
    
    // test: volume
    RealType   Vol_out[4];
    TM.Volume(4, CI_all, Vol_out);
    cout << "Here are the 2-volumes (areas) for these cells: " << endl;
    for (CellIndType ii = 0; ii < 4; ++ii)
        cout << "Cell #" << CI_all[ii] << ",  Volume: " << Vol_out[ii] << endl;
    const bool VOL_CHK = (abs(Vol_out[0]-0.25) > 1E-14) || (abs(Vol_out[1]-0.25) > 1E-14)
                      || (abs(Vol_out[2]-0.25) > 1E-14) || (abs(Vol_out[3]-0.25) > 1E-14);
    if (VOL_CHK)
    {
        cout << "Incorrect!  The volumes (areas) should be 0.25 for all." << endl;
        OUTPUT_CODE = 19;
    }
    cout << endl;
    
    // test: shape regularity
    RealType   SR_out[4];
    TM.Shape_Regularity(4, CI_all, SR_out);
    cout << "Here are the shape regularity ratios for these cells: " << endl;
    for (CellIndType ii = 0; ii < 4; ++ii)
        cout << "Cell #" << CI_all[ii] << ",  SR: " << SR_out[ii] << endl;
    const bool SR_CHK = (abs(SR_out[0]-2.4142135623730950488) > 1E-14) || (abs(SR_out[1]-2.4142135623730950488) > 1E-14)
                     || (abs(SR_out[2]-2.4142135623730950488) > 1E-14) || (abs(SR_out[3]-2.4142135623730950488) > 1E-14);
    if (SR_CHK)
    {
        cout << "Incorrect!  The shape regularity ratios should be ~2.41421356 for all." << endl;
        OUTPUT_CODE = 20;
    }
    cout << endl;
    
    // test: bounding box
    PointType   BB_min[2];
    PointType   BB_max[2];
    TM.Bounding_Box(4, CI_all, BB_min, BB_max);
    cout << "Here is the bounding box of the mesh: " << endl;
    cout << "     BB_min = {" << BB_min[0] << ", " << BB_min[1] << "}," << endl;
    cout << "     BB_max = {" << BB_max[0] << ", " << BB_max[1] << "}." << endl;
    const bool BB_CHK = (abs(BB_min[0]-0.0) > 1E-14) || (abs(BB_min[1]-0.0) > 1E-14)
                     || (abs(BB_max[0]-1.0) > 1E-14) || (abs(BB_max[1]-1.0) > 1E-14);
    if (BB_CHK)
    {
        cout << "Incorrect!  The bounding box should be {0.0, 0.0} to {1.0, 1.0}." << endl;
        OUTPUT_CODE = 21;
    }
    cout << endl;
    
    // test: angles
    RealType   Ang[4*3];
    TM.Angles(4, CI_all, Ang);
    cout << "Here are the (interior) angles, in radians, of each cell in the mesh: " << endl;
    for (CellIndType ii = 0; ii < 4; ++ii)
        cout << "Cell #" << CI_all[ii] << ",  Angles: "
             << Ang[ii*3 + 0] << ", " << Ang[ii*3 + 1] << ", " << Ang[ii*3 + 2] << "." << endl;
    const bool ANG_CHK = (abs(Ang[0]-(PI/2)) > 1E-14) || (abs(Ang[1]-(PI/4)) > 1E-14) || (abs(Ang[2]-(PI/4)) > 1E-14)
                      || (abs(Ang[3]-(PI/2)) > 1E-14) || (abs(Ang[4]-(PI/4)) > 1E-14) || (abs(Ang[5]-(PI/4)) > 1E-14)
                      || (abs(Ang[6]-(PI/2)) > 1E-14) || (abs(Ang[7]-(PI/4)) > 1E-14) || (abs(Ang[8]-(PI/4)) > 1E-14)
                      || (abs(Ang[9]-(PI/2)) > 1E-14) || (abs(Ang[10]-(PI/4)) > 1E-14) || (abs(Ang[11]-(PI/4)) > 1E-14);
    if (ANG_CHK)
    {
        cout << "Incorrect!  The angles should be {pi/2, pi/4, pi/4} for all the cells." << endl;
        OUTPUT_CODE = 22;
    }
    cout << endl;

    
    
    if (OUTPUT_CODE==0)
        cout << "Unit test is successful!" << endl;
    else
    {
        cout << "Unit test failed!" << endl;
        cout << "See OUTPUT_CODE: " << OUTPUT_CODE << endl;
    }
    cout << endl;

    return OUTPUT_CODE;
}

/***/
