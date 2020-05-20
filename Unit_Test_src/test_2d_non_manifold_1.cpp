/*
============================================================================================
   This is a unit test for the AHF data structure implementation.

   This example stores a 2-D non-manifold mesh (consisting of 4 triangles) that has
   1 non-manifold edge and 0 non-manifold vertices.  It also does some basic processing.
   See 'Non_Manifold_Triangle_Mesh_1.jpg' for a picture of the mesh embedded in 3-D.

   Copyright (c) 12-17-2016,  Shawn W. Walker
============================================================================================
*/

#include "../src_code/TypedefMeshes.h"

using namespace std;

// unit test
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, > 0 is failure

    // create the object: manifold 2-D surface mesh example
    SurfaceMesh(MultiMesh);
    // Mesh<3>  MultiMesh(0,0,1,0); // alternative

    // access the parts we need
    BaseMesh<2>& TM = MultiMesh.TriMesh[0];
    BasePtCoord<3>& VX = MultiMesh.Vtx;

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

    // now display coordinates
    cout << endl;
    VX.Display_Vtx_Coord();

    // we now stop adding cells
    TM.Finalize_v2hfs_DEBUG();

    // now display the half-facets (attached to vertices) of the intermediate structure 'v2hfs'
    cout << endl;
    TM.Display_v2hfs();

    // check 'v2hfs' against reference data
    VtxHalfFacetType v2hfs_REF[12];
    v2hfs_REF[ 0].Set(1,0,2);
    v2hfs_REF[ 1].Set(2,0,0);
    v2hfs_REF[ 2].Set(2,0,1);
    v2hfs_REF[ 3].Set(2,1,1);
    v2hfs_REF[ 4].Set(2,2,1);
    v2hfs_REF[ 5].Set(2,3,2);
    v2hfs_REF[ 6].Set(3,1,0);
    v2hfs_REF[ 7].Set(3,1,2);
    v2hfs_REF[ 8].Set(4,2,0);
    v2hfs_REF[ 9].Set(4,2,2);
    v2hfs_REF[10].Set(5,3,0);
    v2hfs_REF[11].Set(5,3,1);
    // error check
    const Vtx2HalfFacet_Mapping& c_V2HF_Map = TM.Get_v2hfs();
    const std::vector<VtxHalfFacetType>& c_v2hfs = c_V2HF_Map.Get_VtxMap();
    for (VtxIndType jj = 0; jj < c_v2hfs.size(); ++jj)
    {
        if (!c_v2hfs[jj].Equal(v2hfs_REF[jj]))
        {
            cout << "Intermediate data 'v2hfs' failed!" << endl;
            OUTPUT_CODE = 1;
            break;
        }
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
    Cell_REF[ 0].Set(1,1,hf.Set());
    Cell_REF[ 0].Set(2,2,hf.Set());
    // Cell #2
    Cell_REF[ 1].Set(0,1,hf.Set());
    Cell_REF[ 1].Set(1,3,hf.Set(2,1));
    Cell_REF[ 1].Set(2,2,hf.Set());
    // Cell #3
    Cell_REF[ 2].Set(0,1,hf.Set());
    Cell_REF[ 2].Set(1,4,hf.Set(3,2));
    Cell_REF[ 2].Set(2,2,hf.Set());
    // Cell #4
    Cell_REF[ 3].Set(0,1,hf.Set());
    Cell_REF[ 3].Set(1,2,hf.Set());
    Cell_REF[ 3].Set(2,5,hf.Set(0,0));
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
    VtxHalfFacetType Vtx2HalfFacets_REF[6];
    Vtx2HalfFacets_REF[ 0].Set( 0,0,2);
    Vtx2HalfFacets_REF[ 1].Set( 1,3,1);
    Vtx2HalfFacets_REF[ 2].Set( 2,3,0);
    Vtx2HalfFacets_REF[ 3].Set( 3,1,2);
    Vtx2HalfFacets_REF[ 4].Set( 4,2,2);
    Vtx2HalfFacets_REF[ 5].Set( 5,3,1);
    // error check
    const std::vector<VtxHalfFacetType>& c_Vtx2HF = TM.Get_Vtx2HalfFacets().Get_VtxMap();
    for (VtxIndType jj = 0; jj < c_Vtx2HF.size(); ++jj)
        if (!c_Vtx2HF[jj].Equal(Vtx2HalfFacets_REF[jj]))
        {
            cout << "'Vtx2HalfFacets' data is incorrect!" << endl;
            OUTPUT_CODE = 3;
            break;
        }

    // display cells attached to a vertex
    VtxIndType V_IN = 2;
    TM.Display_Cells_Attached_To_Vertex(V_IN);
    cout << endl;

    // check attached cells against reference data
    std::vector<CellIndType> cell_ind_1;
    TM.Get_Cells_Attached_To_Vertex(V_IN, 3, cell_ind_1);
    const CellIndType cell_ind_1_REF[4] = {3, 0, 1, 2};
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
    TEST_attached.Set(1,1);
    TM.Display_HalfFacets_Attached_To_HalfFacet(TEST_attached);
    cout << endl;

    // check attached half-facets against reference data
    HalfFacetType attached_REF[4];
    attached_REF[0].Set(1,1);
    attached_REF[1].Set(2,1);
    attached_REF[2].Set(3,2);
    attached_REF[3].Set(0,0);
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
    TEST_attached.Set(3,1);
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
    HalfFacetType non_manifold_hf_REF[1];
    non_manifold_hf_REF[0].Set(3,2);
    std::vector<HalfFacetType> non_manifold_hf;
    TM.Get_Nonmanifold_HalfFacets(non_manifold_hf);
    for (unsigned int jj = 0; jj < non_manifold_hf.size(); ++jj)
        if (!non_manifold_hf[jj].Equal(non_manifold_hf_REF[jj]))
        {
            cout << "Non-manifold facet data is incorrect!" << endl;
            OUTPUT_CODE = 9;
            break;
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
    const VtxIndType uv_REF[6] = {0, 1, 2, 3, 4, 5};
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
    cout << "Vtx #1 and Vtx #5 are connected." << endl;
    if (!TM.Is_Connected(1, 5))
    {
        cout << "Incorrect!  They are NOT connected!" << endl;
        OUTPUT_CODE = 12;
    }
    // test if a pair of vertices is connected by an edge
    cout << "Vtx #0 and Vtx #3 are NOT connected." << endl;
    if (TM.Is_Connected(0, 3))
    {
        cout << "Incorrect!  They ARE connected!" << endl;
        OUTPUT_CODE = 13;
    }
    cout << endl;

    // display all cells attached to an edge
    std::vector<CellIndType> attached_cells;
    MeshEdgeType EE;
    EE.Set(1,2);
    TM.Get_Cells_Attached_To_Edge(EE, attached_cells);
    cout << "Here are the cell indices of cells attached to the edge [1, 2]:" << endl;
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

    if (OUTPUT_CODE==0)
        cout << "Unit test is successful!" << endl;
    else
        cout << "Unit test failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
