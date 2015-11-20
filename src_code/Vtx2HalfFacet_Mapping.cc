/*
============================================================================================
   Class for mapping from a given vertex index to (several) incident half-facets.
   Note: this is generic, meaning this can be used for half-facets in 1-D, 2-D,
         and 3-D meshes.

   EXAMPLE:

      Diagram depicting half-edges (half-facets for 2-D meshes):

                   <2,1>
        V4 +-------------------+ V3
           |\                  |
           |  \          T2    |
           |    \              |
           |      \  <2,2>     |
     <1,2> |        \          | <2,3>
           |    <1,1> \        |
           |            \      |
           |              \    |
           |     T1         \  |
           |                  \|
        V1 +-------------------+ V2
                   <1,3>

   Triangle Connectivity:

   triangle |   vertices
    indices |  V1, V2, V3
   ---------+--------------
       1    |   1,  2,  4
       2    |   2,  3,  4

   Half-Edges attached to vertices:

       Vertex V1:  V1--><1,2>
                   V1--><1,3>
       Vertex V2:  V2--><1,3>
                   V2--><1,1>
                   V2--><2,2>
                   V2--><2,3>
       etc...

   where <ti,ei> is a half-edge attached to V1, where ti and ei define the particular
   half-edge.

   Copyright (c) 11-05-2015,  Shawn W. Walker
============================================================================================
*/

/***************************************************************************************/
/* half-facet data */
struct HalfFacetType
{
    CellIndType   ci; // cell index
    SmallIndType  fi; // local facet (entity) index of ci
    // check equality
    inline bool Equal(const HalfFacetType& IN) const
    {
        const bool ci_CHK = (IN.ci==ci);
        const bool fi_CHK = (IN.fi==fi);
        if (ci_CHK && fi_CHK)
            return true;
        else
            return false;
    }
    // check for NULL
    inline bool Is_Null() const
    {
        if (ci==NULL_IND || fi==NULL_IND)
            return true;
        else
            return false;
    }
    // set value
    inline HalfFacetType& Set(const CellIndType& ci_in, const SmallIndType& fi_in)
    {
        ci = ci_in;
        fi = fi_in;
        return *this;
    }
    inline HalfFacetType& Set(const HalfFacetType& hf_in)
    {
        ci = hf_in.ci;
        fi = hf_in.fi;
        return *this;
    }
    inline HalfFacetType& Set()
    {
        ci = NULL_IND;
        fi = NULL_IND;
        return *this;
    }
};

/***************************************************************************************/
/* vertex/half-facet pair */
struct VtxHalfFacetType
{
    VtxIndType     vtx; // global vertex index
    // half-facet attached to vertex
    CellIndType     ci; // cell index
    SmallIndType    fi; // local facet (entity) index of ci
    // check equality
    inline bool Equal(const VtxHalfFacetType& IN) const
    {
        const bool vtx_CHK = (IN.vtx==vtx);
        const bool  ci_CHK = (IN.ci==ci);
        const bool  fi_CHK = (IN.fi==fi);
        if (vtx_CHK && ci_CHK && fi_CHK)
            return true;
        else
            return false;
    }
    // check for NULL
    inline bool Is_Null() const
    {
        if (vtx==NULL_IND || ci==NULL_IND || fi==NULL_IND)
            return true;
        else
            return false;
    }
    inline VtxHalfFacetType& Set(const VtxIndType& vtx_in, const CellIndType& ci_in, const SmallIndType& fi_in)
    {
        vtx = vtx_in;
        ci  = ci_in;
        fi  = fi_in;
        return *this;
    }
    inline VtxHalfFacetType& Set(const VtxIndType& vtx_in, const HalfFacetType& hf_in)
    {
        vtx = vtx_in;
        ci  = hf_in.ci;
        fi  = hf_in.fi;
        return *this;
    }
    inline VtxHalfFacetType& Set()
    {
        vtx = NULL_IND;
        ci  = NULL_IND;
        fi  = NULL_IND;
        return *this;
    }
	// inline bool operator == (const VtxHalfFacetType& vhf) const
    // {
        // return Equal(vhf);
    // }
};

/* C++ class definition */
#define  V2HF  Vtx2HalfFacet_Mapping
class V2HF
{
public:
    V2HF();
    ~V2HF();
    void Clear() { VtxMap.clear(); };
    // get number of entries in the Vertex-to-HalfFacet map
    VtxIndType Size() const { return (VtxIndType) VtxMap.size(); };
    // return const reference to internal data
    const std::vector<VtxHalfFacetType>& Get_VtxMap() const { return VtxMap; };
    // allocate space in the Vertex-to-HalfFacet map
    void Reserve(const VtxIndType&);

    // append <vertex,half-facet> pair
    void Append(const VtxIndType&, const HalfFacetType&);
    void Append(const VtxIndType&, VtxHalfFacetType&);
    void Append(const VtxHalfFacetType&);
    // sort the VtxMap so it is useable
    void Sort();

    /* Note: all public methods below this line require VtxMap to be sorted
       before they will work correctly.  Make sure to run 'Sort()' first! */
    // find one half-facet attached to the given vertex
    unsigned int Get_Half_Facet(const VtxIndType&, std::vector<VtxHalfFacetType>::const_iterator&) const;
    unsigned int Get_Half_Facet(const VtxIndType&, HalfFacetType&) const;
    // find all half-facets
    unsigned int Get_Half_Facets(const VtxIndType&, std::pair <std::vector<VtxHalfFacetType>::const_iterator,
                                                               std::vector<VtxHalfFacetType>::const_iterator>&) const;
    unsigned int Get_Half_Facets(const VtxIndType&, std::pair <std::vector<VtxHalfFacetType>::iterator,
                                                               std::vector<VtxHalfFacetType>::iterator>&);
    unsigned int Get_Half_Facets(const VtxIndType&, std::vector<HalfFacetType>&) const;

    // print out the half-facets attached to given vertex
    void Display_Half_Facets(const VtxIndType&) const;
    // get a vector of *unique* vertices from VtxMap
    void Get_Unique_Vertices(std::vector<VtxIndType>&) const;
    void Display_Unique_Vertices() const;

private:
    // map from vertices to adjacent half-facets
    /* note: the same vertex can appear multiple times, but with different attached half-facets. */
    std::vector<VtxHalfFacetType> VtxMap; // with duplicate vertices!

    // number between 0.0 and 1.0 denoting percent extra space to allocate
    double  Reserve_Buffer;
};

/***************************************************************************************/
/* constructor */
V2HF::V2HF()
{
    // ensure memory is clear to start
    Clear();
    Reserve_Buffer = 0.2;
}

/***************************************************************************************/
/* DE-structor */
V2HF::~V2HF()
{
    // clear the data
    Clear();
}

/***************************************************************************************/
/* Allocate memory to hold a bunch of vertex-to-half-facet structs of given size
  (plus a little). */
void V2HF::Reserve(const VtxIndType& Num_V2HF)
{
    // compute the desired size (with extra) to allocate
    VtxIndType Desired_SIZE = (VtxIndType) ((1.0 + Reserve_Buffer) * Num_V2HF);
    VtxMap.reserve(Desired_SIZE);
}

/***************************************************************************************/
/* append a <vertex, half-facet> pair. */
void V2HF::Append(const VtxIndType& vi, const HalfFacetType& hf)
{
    VtxHalfFacetType TEMP;
    TEMP.vtx = vi;
    TEMP.ci  = hf.ci;
    TEMP.fi  = hf.fi;
    Append(TEMP);
}
void V2HF::Append(const VtxIndType& vi, VtxHalfFacetType& vhf)
{
    vhf.vtx = vi; // the rest of the struct is ready!
    Append(vhf);
}
void V2HF::Append(const VtxHalfFacetType& in)
{
    const VtxIndType current_size = (VtxIndType) VtxMap.size();
    if (current_size == VtxMap.capacity())
        Reserve(current_size);
    VtxMap.push_back(in);
}

/***************************************************************************************/
/* sort the <vertex, half-facet> list. */
bool VtxMap_order (const VtxHalfFacetType& Va, const VtxHalfFacetType& Vb)
{
	if (Va.vtx==Vb.vtx)
	{
		if (Va.ci==Vb.ci)
			return (Va.fi < Vb.fi); 
		else
			return (Va.ci < Vb.ci);
	}
	else
		return (Va.vtx < Vb.vtx);
}
// (SWW) Note: code can be made ~%7 faster if VtxMap_order is simplified to "return (Va.vtx < Vb.vtx);"
void V2HF::Sort()
{
    std::sort(VtxMap.begin(), VtxMap.end(), VtxMap_order);
}

/***************************************************************************************/
/* find a single half-facet attached to the given vertex. */
unsigned int V2HF::Get_Half_Facet(const VtxIndType& vi, std::vector<VtxHalfFacetType>::const_iterator& it) const
{
    std::pair <std::vector<VtxHalfFacetType>::const_iterator,
               std::vector<VtxHalfFacetType>::const_iterator> RR;
    const unsigned int Num_HF = Get_Half_Facets(vi, RR);

    if (Num_HF==0) // nothing was found, so have iterator point to the end
        it = VtxMap.end();
    else // at least one was found, so point to the first
        it = RR.first;

    return Num_HF;
}
/* find a single half-facet attached to the given vertex, and return half-facet. */
unsigned int V2HF::Get_Half_Facet(const VtxIndType& vi, HalfFacetType& hf) const
{
    std::vector<VtxHalfFacetType>::const_iterator it;
    const unsigned int Num_HF = Get_Half_Facet(vi, it);

    if (it!=VtxMap.end()) // found valid key
    {
        // copy over
        hf.ci = (*it).ci;
        hf.fi = (*it).fi;
    }
    else
    {
        // set to null value
        hf.Set();
    }
    return Num_HF;
}

/***************************************************************************************/
/* find *all* half-facets attached to the given vertex. */
bool VtxMap_vtx_equal(const VtxHalfFacetType& Va, const VtxHalfFacetType& Vb) { return (Va.vtx < Vb.vtx); }
// note: consider the "VtxHalfFacetType"'s to be equal if the vertex indices match.
unsigned int V2HF::Get_Half_Facets(const VtxIndType& vi, std::pair <std::vector<VtxHalfFacetType>::const_iterator,
                                                                    std::vector<VtxHalfFacetType>::const_iterator>& range) const
{
    VtxHalfFacetType  TEMP;
    TEMP.vtx = vi;
    range = std::equal_range(VtxMap.begin(), VtxMap.end(), TEMP, VtxMap_vtx_equal);
    const unsigned int Num_HF = (unsigned int) std::distance(range.first,range.second);
    return Num_HF;
}
/* non-const version of above method! */
unsigned int V2HF::Get_Half_Facets(const VtxIndType& vi, std::pair <std::vector<VtxHalfFacetType>::iterator,
                                                                    std::vector<VtxHalfFacetType>::iterator>& range)
{
    VtxHalfFacetType  TEMP;
    TEMP.vtx = vi;
    range = std::equal_range(VtxMap.begin(), VtxMap.end(), TEMP, VtxMap_vtx_equal);
    const unsigned int Num_HF = (unsigned int) std::distance(range.first,range.second);
    return Num_HF;
}
/* find *all* half-facets attached to the given vertex, and return in a std::vector.
   Note: in most cases, this routine is not needed (more for testing). */
unsigned int V2HF::Get_Half_Facets(const VtxIndType& vi, std::vector<HalfFacetType>& hf) const
{
    std::pair <std::vector<VtxHalfFacetType>::const_iterator,
               std::vector<VtxHalfFacetType>::const_iterator>  range;
    const unsigned int Num_HF = Get_Half_Facets(vi, range);

    hf.clear(); // start clean!
    if (Num_HF > 0) // at least one half-facet was found
    {
        hf.reserve(Num_HF);
        for (std::vector<VtxHalfFacetType>::const_iterator it=range.first; it!=range.second; ++it)
        {
            // append
            HalfFacetType  TEMP;
            TEMP.ci = (*it).ci;
            TEMP.fi = (*it).fi;
            hf.push_back(TEMP);
        }
    }
    // else // key not found
    return Num_HF;
}

/***************************************************************************************/
/* display half-facets attached to a given vertex.
   Note: you must run Sort() before using this! */
void V2HF::Display_Half_Facets(const VtxIndType& vi=NULL_IND) const
{
    if (vi==NULL_IND) // assume we want all the half-facet(s) for all the stored vertices
    {
        std::vector<VtxHalfFacetType>::const_iterator it = VtxMap.begin();
        VtxIndType Prev_Vtx, Current_Vtx;
        Prev_Vtx = NULL_IND; // init
        std::cout << "Vertex and attached half-facets, <cell index, local facet index>:";
        while (it!=VtxMap.end())
        {
            // get the current vertex
            Current_Vtx = (*it).vtx;
            if (Current_Vtx!=Prev_Vtx)
            {
                // found a new vertex
                std::cout << std::endl;
                std::cout << "Vtx# " << Current_Vtx << ": ";
                std::cout << "<" << (*it).ci << ", " <<  (*it).fi << ">";
                Prev_Vtx = Current_Vtx; // update
            }
            else if (Current_Vtx==Prev_Vtx) // still with the same vertex
            {
                // print the current half-facet
                std::cout << ", <" << (*it).ci << ", " <<  (*it).fi << ">";
            }
            ++it; // go to the next half-facet
        }
    }
    else // print out half-facets for one particular vertex
    {
        std::vector<HalfFacetType> HF1;
        Get_Half_Facets(vi, HF1);
        std::cout << std::endl;
        std::cout << "Half-facets <cell index, local facet index> attached to Vtx# " << vi << ":" << std::endl;
        for (VtxIndType kk = 0; kk < HF1.size(); kk++)
            std::cout << "<" << HF1[kk].ci << ", " <<  HF1[kk].fi << ">, ";
    }
    std::cout << std::endl;
}

/***************************************************************************************/
/* get unique list of vertices (VtxMap should already be sorted).
   Note: this just prints info out. */
void V2HF::Get_Unique_Vertices(std::vector<VtxIndType>& unique_vertices) const
{
    // extract all the vertex indices
    unique_vertices.clear(); // start fresh
    unique_vertices.reserve(Size());
    for (std::vector<VtxHalfFacetType>::const_iterator it=VtxMap.begin(); it!=VtxMap.end(); ++it)
        unique_vertices.push_back((*it).vtx);

    // get unique set of vertices
    std::vector<VtxIndType>::iterator it_end;
    it_end = std::unique_copy(unique_vertices.begin(), unique_vertices.end(), unique_vertices.begin());
    // get number of unique vertices
    const VtxIndType LENGTH = (unsigned int) std::distance(unique_vertices.begin(),it_end);
    // resize
    unique_vertices.resize(LENGTH);
}
/* just print them to the screen. */
void V2HF::Display_Unique_Vertices() const
{
    // extract all the vertex indices
    std::vector<VtxIndType> unique_vertices;
    Get_Unique_Vertices(unique_vertices);

    // now print them out
    std::cout << "Unique list of vertex indices: " << std::endl;
    std::cout << *(unique_vertices.begin());
    for (std::vector<VtxIndType>::const_iterator vi=unique_vertices.begin()+1; vi!=unique_vertices.end(); ++vi)
    {
        std::cout << ", " << (*vi);
    }
    std::cout << std::endl;
}

#undef V2HF

/***/
