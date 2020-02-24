/*
============================================================================================
   This is a demo for how to use the AHF Mesh class.

   We create a mesh of 5 tetrahedra, with one facet having 3 tetrahedra attached.
   And an edge of the mesh has only two tetrahedra attached (in a non-manifold way).
   See the .jpg file in this directory to see a picture.

   The connectivity data is given by:
   
    cells |  vertices  |  neighboring half-facets
      C0  | 0, 1, 2, 3 | <1,3>, <-,->, <-,->, <-,->
      C1  | 1, 2, 3, 4 | <-,->, <-,->, <-,->, <2,1>
      C2  | 1, 5, 2, 3 | <-,->, <0,0>, <-,->, <-,->
      C3  | 1, 8, 6, 5 | <4,2>, <-,->, <-,->, <-,->
      C4  | 6, 8, 7, 5 | <-,->, <-,->, <3,0>, <-,->

   where <Ci,Fi> denotes the cell index and the local facet index of that cell.
   
   This example stores the mesh and does some basic processing.

   Note: see "BaseMesh.cc" comments for more info.

   Copyright (c) 02-11-2020,  Shawn W. Walker
============================================================================================
*/

#include "../../src_code/TypedefMeshes.h"

using namespace std;

// demo
int main()
{
    // init output code
    int OUTPUT_CODE = 0; // 0 indicates success, > 0 is failure

    // create the object: a 2-D topological mesh
    BaseMesh<3>  BM;

    // define the cell connectivity (5 cells)
    BM.Reserve_Cells(5);
    BM.Append_Cell(0,1,2,3);
    BM.Append_Cell(1,2,3,4);
    BM.Append_Cell(1,5,2,3);
    BM.Append_Cell(1,8,6,5);
    BM.Append_Cell(6,8,7,5);
	cout << endl;

	// now build internal connectivity information
	BM.Finalize_Mesh_Connectivity();

    // let's not risk changing the mesh inadvertently
    const BaseMesh<3>& BM_c = BM;

	// and display cell connectivity data
    cout << endl;
    BM_c.Display_Cell();
    cout << endl;

    // extract cell #k
    const CellIndType k0 = 2;
    const CellSimplexType<3>& Ck = BM_c.Get_Cell(k0);
    const VtxIndType* Ck_vtx = Ck.Get_vtx();
    const HalfFacetType* Ck_hf = Ck.Get_halffacet();
    // and here is the info:
	cout << "The global vertices of Cell #" << k0 << ": ";
    cout << Ck_vtx[0] << ", " << Ck_vtx[1] << ", " << Ck_vtx[2] << ", " << Ck_vtx[3] << endl;
	cout << "Cell #" << k0 << " has four local facets: F0, F1, F2, F3" << endl;
    cout << "      F0 is identical with local face [V1, V2, V3]: [" <<
                   Ck_vtx[1] << ", " << Ck_vtx[2] << ", " << Ck_vtx[3] << "]" << endl;
    cout << "      F1 is identical with local face [V0, V3, V2]: [" <<
                   Ck_vtx[0] << ", " << Ck_vtx[3] << ", " << Ck_vtx[2] << "]" << endl;
    cout << "      F2 is identical with local face [V0, V1, V3]: [" <<
                   Ck_vtx[0] << ", " << Ck_vtx[1] << ", " << Ck_vtx[3] << "]" << endl;
    cout << "      F3 is identical with local face [V0, V2, V1]: [" <<
                   Ck_vtx[0] << ", " << Ck_vtx[2] << ", " << Ck_vtx[1] << "]" << endl;
	cout << "The neighboring half-facet of F0 is ..." << endl;
    if (Ck_hf[0].Is_Null())
    {
        cout << "      a NULL half-facet." << endl;
    }
    else
    {
        cout << "      contained in Cell index: " << Ck_hf[0].ci << ", with local facet index: " << Ck_hf[0].fi << endl;
        cout << "      But this is incorrect!" << endl;
        OUTPUT_CODE = 1;
    }
	cout << "The neighboring half-facet of F1 is ..." << endl;
    cout << "      contained in Cell index: " << Ck_hf[1].ci << ", with local facet index: " << Ck_hf[1].fi << endl;
    cout << "The neighboring half-facet of F2 is ..." << endl;
    if (Ck_hf[2].Is_Null())
    {
        cout << "      a NULL half-facet." << endl;
    }
    else
    {
        cout << "      contained in Cell index: " << Ck_hf[2].ci << ", with local facet index: " << Ck_hf[2].fi << endl;
        cout << "      But this is incorrect!" << endl;
        OUTPUT_CODE = 2;
    }
    cout << "The neighboring half-facet of F3 is ..." << endl;
    if (Ck_hf[3].Is_Null())
    {
        cout << "      a NULL half-facet." << endl;
    }
    else
    {
        cout << "      contained in Cell index: " << Ck_hf[3].ci << ", with local facet index: " << Ck_hf[3].fi << endl;
        cout << "      But this is incorrect!" << endl;
        OUTPUT_CODE = 3;
    }
    cout << endl;

    // can also view the vertex-to-halffacet mapping
    // NOTE: facets are faces for a 3-D topological mesh
    BM_c.Display_Vtx2HalfFacets();
    cout << endl;

    // can also get the vertex-to-halffacet mapping explicitly
    const Vtx2HalfFacet_Mapping V2HF = BM_c.Get_Vtx2HalfFacets();
    const std::vector<VtxHalfFacetType> vtx_hf_vec = V2HF.Get_VtxMap();
	// output the data "manually"
	cout << "The vertex-to-halffacet mapping is (manual):" << endl;
    for (std::vector<VtxHalfFacetType>::const_iterator it=vtx_hf_vec.begin(); it!=vtx_hf_vec.end(); ++it)
    {
        cout << "Vtx #" << (*it).vtx << " is attached to half-facet: ";
        cout << "<" << (*it).ci << ", " << (*it).fi << ">" << endl;
    }
	cout << endl;

    // display the cells attached to vertex #1
    const VtxIndType v1 = 1;
    BM_c.Display_Cells_Attached_To_Vertex(v1);
    cout << endl;

    // can also get a vector of the attached cell indices
	std::vector<CellIndType> cell_ind_1;
	BM_c.Get_Cells_Attached_To_Vertex(v1, cell_ind_1);
	// output the cell indices "manually"
	cout << "The attached cells to vertex # " << v1 << " are (manual):" << endl;
    for (std::vector<CellIndType>::const_iterator it=cell_ind_1.begin(); it!=cell_ind_1.end(); ++it)
    {
        cout << *it << endl;
    }
	cout << endl;

    // display half-facets attached to a given half-facet
    HalfFacetType given_hf;
    given_hf.Set(2,1); // half-facet is <cell #2, local facet #1>
    BM_c.Display_HalfFacets_Attached_To_HalfFacet(given_hf);
    cout << endl;

    // can also get a vector of the attached half-facets
    std::vector<HalfFacetType> attached;
    BM_c.Get_HalfFacets_Attached_To_HalfFacet(given_hf, attached);
	// output the half-facets "manually"
	cout << "The half-facets attached to half-facet (manual): ";
	given_hf.Print();
	cout << " are: " << endl;
    for (std::vector<HalfFacetType>::const_iterator it=attached.begin(); it!=attached.end(); ++it)
	{
        (*it).Print();
		cout << endl;
	}
	cout << endl;

    // check that half-facet has no neighbors
    given_hf.Set(3,3); // half-facet is <cell #3, local facet #3>
    BM_c.Get_HalfFacets_Attached_To_HalfFacet(given_hf, attached);
    if (attached.size()!=1)
    {
        cout << "HalfFacet attachment data is incorrect!" << endl;
        OUTPUT_CODE = 4;
    }
    else
    {
        cout << "This half-facet: " << endl;
        given_hf.Print();
        cout << "    has no neighbor." << endl;
    }
	cout << endl;

    // display all edges in the mesh
    BM_c.Display_Edges();

    // can also get a vector of the edges
    std::vector<MeshEdgeType> edges;
    BM_c.Get_Edges(edges);
    // output the edges "manually"
    cout << "A unique list of edges in the mesh (manual):" << endl;
    for (std::vector<MeshEdgeType>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
    {
        cout << "[" << (*it).vtx[0] << ", " << (*it).vtx[1] << "]" << endl;
    }
    cout << endl;

    // test if a pair of vertices is connected by an edge
    VtxIndType vv[5] = {0, 2, 6, 4, 3};
    for (VtxIndType ii = 0; ii < 4; ++ii)
    {
        cout << "Vtx #" << vv[ii] << " and Vtx #" << vv[ii+1] << " are ";
        if (BM_c.Is_Connected(vv[ii], vv[ii+1]))
            cout << "connected." << endl;
        else
            cout << "NOT connected." << endl;
    }
    cout << endl;

    // display all cells attached to an edge
    std::vector<CellIndType> attached_cells;
    MeshEdgeType EE;
    EE.Set(1,5); // edge defined by [V1, V5] using global vertex indices
    BM_c.Get_Cells_Attached_To_Edge(EE, attached_cells);
    cout << "Here are the cell indices of cells attached to the edge ["
             << EE.vtx[0] << ", " << EE.vtx[1] << "]:" << endl;
    std::vector<CellIndType>::const_iterator ic;
    for (ic = attached_cells.begin(); ic!=attached_cells.end(); ++ic)
        cout << "Cell #" << (*ic) << endl;
    cout << endl;

    // get the free boundary
	// i.e. get the set of half-facets that are attached to only *one* cell
    std::vector<HalfFacetType> free_bdy;
    BM_c.Get_FreeBoundary(free_bdy);
    cout << "Here are the half-facets that lie on the free boundary:" << endl;
    for (std::vector<HalfFacetType>::const_iterator hfi = free_bdy.begin(); hfi!=free_bdy.end(); ++hfi)
    {
        (*hfi).Print();
        cout << endl;
    }
    cout << endl;

    // display non-manifold half-facets
    BM_c.Display_Nonmanifold_HalfFacets();
    cout << endl;

    // can also get a vector of the non-manifold half-facets
    std::vector<HalfFacetType> non_manifold_hf;
    BM_c.Get_Nonmanifold_HalfFacets(non_manifold_hf);
	// output the half-facets "manually"
	cout << "The non-manifold half-facets are (manual): ";
    for (std::vector<HalfFacetType>::const_iterator it=non_manifold_hf.begin(); it!=non_manifold_hf.end(); ++it)
	{
        (*it).Print();
		cout << endl;
	}
	cout << endl;
    if (non_manifold_hf.size() > 1)
    {
        cout << "There should only be one non-manifold half-facet!" << endl;
        OUTPUT_CODE = 5;
    }
    /* one can now call "Get_HalfFacets_Attached_To_HalfFacet"
           to see what all comprises the non-manifold half-facet. */

    // display non-manifold vertices
    BM_c.Display_Nonmanifold_Vertices();
    cout << endl;

    // can also get a vector of the non-manifold vertices
	std::vector<VtxIndType> non_man_vtx;
    BM_c.Get_Nonmanifold_Vertices(non_man_vtx);
	// output the vertex indices "manually"
	cout << "The non-manifold vertices are (manual):" << endl;
    for (std::vector<VtxIndType>::const_iterator it=non_man_vtx.begin(); it!=non_man_vtx.end(); ++it)
    {
        cout << *it << endl;
    }
	cout << endl;
    if (non_man_vtx.size()!=2)
    {
        cout << "There should only be two non-manifold vertices!" << endl;
        OUTPUT_CODE = 6;
    }

    if (OUTPUT_CODE==0)
        cout << "Demo completed successfully!" << endl;
    else
        cout << "Demo failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
