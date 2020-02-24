/*
============================================================================================
   This is a demo for how to use the AHF Mesh class.

   We create a mesh of 5 triangles, with one edge having 4 triangles attached
   in a "fan" pattern.  Embedding it into \R^3, it looks like:

                                        +V3
                                       /|
                                      / |
                                     /  |
                               <-,->/   |
                                   /    |<-,->
                                  /     |
                V4    <-,->    V1/<1,0> |     <-,->
                 +--------------+       |-------------+V2
                  \__   <2,0>    \ <1,2>|     <0,0>  /
                     \__          \     | C0       /
                 <-,->  \__<2,1>   \ C1 |   <0,2>/
                           \__   C2 \   |      / <-,->
                           /  \__    \  |    /
                          /      \__  \ |  /
                         /          \__\|/
                        /<3,0>          +V0
                       /             __/ \
                 <-,->/      <3,1>__/     \
                     /  C3    ___/    <4,1>\<-,->
                    /     ___/<4,2>         \
                   /  ___/           C4      \
                  /__/        <4,0>           \
                 +-----------------------------+
                V5            <-,->           V6

   where the Vi are the vertex indices, and Ci are the cell indices;
   some of the half-facets are also shown.  Obviously, this is not a
   manifold mesh.  The half-facets attached to the edge [V0, V1]
   (not shown in diagram) are:  <0,1>, <1,1>, <2,2>, <3,2>.

   This example stores the mesh and does some basic processing.

   Note: see "BaseMesh.cc" comments for more info.

   Copyright (c) 02-03-2020,  Shawn W. Walker
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
    BaseMesh<2>  BM;

    // define the cell connectivity (5 cells)
    BM.Reserve_Cells(5);
    BM.Append_Cell(0,2,1);
    BM.Append_Cell(0,3,1);
    BM.Append_Cell(0,1,4);
    BM.Append_Cell(0,1,5);
    BM.Append_Cell(0,5,6);
	cout << endl;

	// now build internal connectivity information
	BM.Finalize_Mesh_Connectivity();

    // let's not risk changing the mesh inadvertently
    const BaseMesh<2>& BM_c = BM;

	// and display cell connectivity data
    cout << endl;
    BM_c.Display_Cell();
    cout << endl;

    // extract cell #k
    const CellIndType k0 = 3;
    const CellSimplexType<2>& Ck = BM_c.Get_Cell(k0);
    const VtxIndType* Ck_vtx = Ck.Get_vtx();
    const HalfFacetType* Ck_hf = Ck.Get_halffacet();
    // and here is the info:
	cout << "The global vertices of Cell #" << k0 << ": ";
    cout << Ck_vtx[0] << ", " << Ck_vtx[1] << ", " << Ck_vtx[2] << endl;
	cout << "Cell #" << k0 << " has three local facets: F0, F1, F2" << endl;
    cout << "      F0 is identical with local edge [V1, V2]: [" << Ck_vtx[1] << ", " << Ck_vtx[2] << "]" << endl;
    cout << "      F1 is identical with local edge [V2, V0]: [" << Ck_vtx[2] << ", " << Ck_vtx[0] << "]" << endl;
    cout << "      F2 is identical with local edge [V0, V1]: [" << Ck_vtx[0] << ", " << Ck_vtx[1] << "]" << endl;
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
    cout << "      contained in Cell index: " << Ck_hf[2].ci << ", with local facet index: " << Ck_hf[2].fi << endl;
    cout << endl;

    // can also view the vertex-to-halffacet mapping
    // NOTE: facets are edges for a 2-D topological mesh
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

    // display the cells attached to vertex #0
    const VtxIndType v0 = 0;
    BM_c.Display_Cells_Attached_To_Vertex(v0);
    cout << endl;

    // can also get a vector of the attached cell indices
	std::vector<CellIndType> cell_ind_1;
	BM_c.Get_Cells_Attached_To_Vertex(v0, cell_ind_1);
	// output the cell indices "manually"
	cout << "The attached cells to vertex # " << v0 << " are (manual):" << endl;
    for (std::vector<CellIndType>::const_iterator it=cell_ind_1.begin(); it!=cell_ind_1.end(); ++it)
    {
        cout << *it << endl;
    }
	cout << endl;

    // display half-facets attached to a given half-facet
    HalfFacetType given_hf;
    given_hf.Set(2,2); // half-facet is <cell #2, local facet #2>
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
    given_hf.Set(4,1); // half-facet is <cell #4, local facet #1>
    BM_c.Get_HalfFacets_Attached_To_HalfFacet(given_hf, attached);
    if (attached.size()!=1)
    {
        cout << "HalfFacet attachment data is incorrect!" << endl;
        OUTPUT_CODE = 2;
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
    VtxIndType vv[4] = {0, 3, 4, 6};
    for (VtxIndType ii = 0; ii < 3; ++ii)
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
    EE.Set(0,5); // edge defined by [V0, V5] using global vertex indices
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
        OUTPUT_CODE = 3;
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
    if (non_man_vtx.size() > 0)
    {
        cout << "There should be NO non-manifold vertices!" << endl;
        OUTPUT_CODE = 4;
    }

    if (OUTPUT_CODE==0)
        cout << "Demo completed successfully!" << endl;
    else
        cout << "Demo failed!" << endl;
    cout << endl;

    return OUTPUT_CODE;
}

/***/
