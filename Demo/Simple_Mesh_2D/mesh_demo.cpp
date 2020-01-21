/*
============================================================================================
   This is a demo for how to use the AHF Mesh class.

   This example stores a 2-D manifold mesh (consisting of 4 triangles).
   It also does some basic processing.

   Copyright (c) 12-19-2016,  Shawn W. Walker
============================================================================================
*/

#include "../../src_code/TypedefMeshes.h"

using namespace std;

// unit test
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, 1 is failure

    // create the object: a 2-D mesh in \R^2 (the x-y plane)
    TriMesh(MultiMesh);
    // Mesh<2>  MultiMesh(0,1,0); // alternative
    
    // access the parts we need
    BaseMesh<2>& TM = MultiMesh.TriMesh[0];
    BasePtCoord<2>& VX = MultiMesh.Vtx;

    // define the cell connectivity (4 cells)
    TM.Reserve_Cells(4);
    TM.Append_Cell(0,1,4);
    TM.Append_Cell(1,2,4);
    TM.Append_Cell(2,3,4);
    TM.Append_Cell(3,0,4);

    // now add the vertex point coordinates (5 vertices)
    VX.Init_Points(5);
    VX.Set_Coord(0, 0.0,0.0);
    VX.Set_Coord(1, 1.0,0.0);
    VX.Set_Coord(2, 1.0,1.0);
    VX.Set_Coord(3, 0.0,1.0);
    VX.Set_Coord(4, 0.5,0.5);

    // display coordinates
    VX.Display_Vtx_Coord();

	// note: mesh is a square composed of 4 triangles
	
	// now build internal connectivity information
    cout << endl;
	TM.Finalize_Mesh_Connectivity();
	
	// and display cell connectivity data
    cout << endl;
    TM.Display_Cell();
    cout << endl;

    // can also view the vertex-to-halffacet mapping
    TM.Display_Vtx2HalfFacets();
    cout << endl;

    // display the cells attached to vertex #2
    VtxIndType V_IN = 2;
    TM.Display_Cells_Attached_To_Vertex(V_IN);
    cout << endl;

    // can also get a vector of the attached cell indices
	std::vector<CellIndType> cell_ind_1;
	TM.Get_Cells_Attached_To_Vertex(V_IN, cell_ind_1);
	// output the cell indices "manually"
	cout << "The attached cells to vertex # " << V_IN << " are:" << endl;
    for (std::vector<CellIndType>::const_iterator it=cell_ind_1.begin(); it!=cell_ind_1.end(); ++it)
        cout << "#" << *it << endl;
	cout << endl;

    // display half-facets attached to a given half-facet
    HalfFacetType given_hf;
    given_hf.Set(3,1); // half-facet is <cell #3, local facet #1>
    TM.Display_HalfFacets_Attached_To_HalfFacet(given_hf);
    cout << endl;

    // can also get a vector of the attached half-facets
    std::vector<HalfFacetType> attached;
    TM.Get_HalfFacets_Attached_To_HalfFacet(given_hf, attached);
	// output the half-facets "manually"
	cout << "The half-facets attached to half-facet: ";
	given_hf.Print();
	cout << " are: " << endl;
    for (std::vector<HalfFacetType>::const_iterator it=attached.begin(); it!=attached.end(); ++it)
	{
        (*it).Print();
		cout << endl;
	}
	cout << endl;

    // check that half-facet has no neighbors
    given_hf.Set(2,2); // half-facet is <cell #2, local facet #2>
    TM.Get_HalfFacets_Attached_To_HalfFacet(given_hf, attached);
    if (attached.size()!=1)
    {
        cout << "HalfFacet attachment data is incorrect!" << endl;
        OUTPUT_CODE = 1;
    }
	
    // display all edges in the mesh
    std::vector<MeshEdgeType> edges;
    TM.Get_Edges(edges);
    cout << "A unique list of edges in the mesh:" << endl;
    for (std::vector<MeshEdgeType>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
    {
        cout << "[" << (*it).vtx[0] << ", " << (*it).vtx[1] << "]" << endl;
    }
    cout << endl;

    // test if a pair of vertices is connected by an edge
    cout << "Vtx #0 and Vtx #4 are connected." << endl;
    if (!TM.Is_Connected(0, 4))
    {
        cout << "Incorrect!  They are NOT connected!" << endl;
        OUTPUT_CODE = 2;
    }
    // test if a pair of vertices is connected by an edge
    cout << "Vtx #0 and Vtx #2 are NOT connected." << endl;
    if (TM.Is_Connected(0, 2))
    {
        cout << "Incorrect!  They ARE connected!" << endl;
        OUTPUT_CODE = 3;
    }
    cout << endl;

    // display all cells attached to an edge
    std::vector<CellIndType> attached_cells;
    MeshEdgeType EE;
    EE.Set(1,4); // edge defined by [V1, V4] using global vertex indices
    TM.Get_Cells_Attached_To_Edge(EE, attached_cells);
    cout << "Here are the cell indices of cells attached to the edge [1, 4]:" << endl;
    std::vector<CellIndType>::const_iterator ic;
    for (ic = attached_cells.begin(); ic!=attached_cells.end(); ++ic)
        cout << "Cell #" << (*ic) << endl;
    cout << endl;

    // get the free boundary
	// i.e. get the set of half-facets that are attached to only *one* cell
    std::vector<HalfFacetType> free_bdy;
    TM.Get_FreeBoundary(free_bdy);
    cout << "Here are the half-facets that lie on the free boundary:" << endl;
    for (std::vector<HalfFacetType>::const_iterator hfi = free_bdy.begin(); hfi!=free_bdy.end(); ++hfi)
    {
        (*hfi).Print();
        cout << endl;
    }
    cout << endl;

    // display non-manifold facets
    TM.Display_Nonmanifold_HalfFacets();
    cout << endl;

    // display non-manifold vertices
    TM.Display_Nonmanifold_Vertices();
    cout << endl;

    if (OUTPUT_CODE==0)
        cout << "Demo completed successfully!" << endl;
    else
        cout << "Demo failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
