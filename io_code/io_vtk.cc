/*
============================================================================================
   Class for reading/writing mesh data to a .vtk file format.  Can do both ASCII and
   BINARY format.
   
   Note: this routine is limited to simplex meshes: line-segment, triangle, tetrahedron.

   Copyright (c) 12-15-2016,  Shawn W. Walker
============================================================================================
*/

#define _IO_VTK_CC

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

/* C++ class definition */
#define  MESH_IO  Mesh_IO_VTK
class MESH_IO
{
public:
    MESH_IO();
    ~MESH_IO();
	
	// write the mesh to a VTK file
	template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
	void Write_ASCII_VTK(std::string, std::string, const Mesh<CELL_DIM, GEO_DIM>&);
	template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
	void Write_Binary_VTK(std::string, std::string, const Mesh<CELL_DIM, GEO_DIM>&);
	
	// read a VTK file to determine the type of mesh and store RAW data
	MeshInterface* Read_ASCII_VTK(const std::string&);
	MeshInterface* Read_BINARY_VTK(const std::string&);

	// convert raw data to a Mesh object (using MeshInterface)
	MeshInterface* Convert_to_Mesh_Object(const SmallIndType&, const SmallIndType&,
									      const CellIndType&, const VtxIndType&,
									      const std::vector<CellIndType>&,
										  const std::vector<PointType>&);

private:

    // display raw mesh stats
    void Display_Mesh_Info(const std::string&, const std::string&, const SmallIndType&, const SmallIndType&,
                           const CellIndType&, const VtxIndType&);

	// simple write routines
	void Write_VTK_Hdr(std::ofstream&, const std::string&, const std::string&);
	void Write_VTK_Vtx_Hdr(std::ofstream&, const VtxIndType&, const std::string&);
	SmallIndType Write_VTK_Cell_Data_Hdr (std::ofstream&, const SmallIndType&, const CellIndType&);
    SmallIndType Write_VTK_Cell_Types_Hdr(std::ofstream&, const SmallIndType&, const CellIndType&);
	
	// simple read routines
	VtxIndType Read_VTK_Initial(std::ifstream&, const std::string&, const bool&);
	void Read_VTK_Cell_Data_Hdr(std::ifstream&, const std::string&, SmallIndType&, CellIndType&);
	void Read_VTK_Cell_Type_Hdr(std::ifstream&, const std::string&, const CellIndType&);
	
	// helper
	SmallIndType Cell_Label_to_Top_Dim(const SmallIndType&, std::string&) const;
	
	// conversion routines to deal with Big vs. Little Endian stuff
	bool System_Is_Little_Endian();
	void PointTypeFl_to_LitEnd_Char(const PointTypeFl& I, char O[4]);
	void PointTypeFl_to_BigEnd_Char(const PointTypeFl& I, char O[4]);
	void PointType_to_LitEnd_Char(const PointType& I, char O[8]);
	void PointType_to_BigEnd_Char(const PointType& I, char O[8]);
	void BigEnd_Char_to_LitEnd_Char_32bit(const char I[4], char O[4]);
	void BigEnd_Char_to_LitEnd_Char_64bit(const char I[8], char O[8]);
	
	void VtxIndType_to_LitEnd_Char(const VtxIndType& I, char O[4]);
	void VtxIndType_to_BigEnd_Char(const VtxIndType& I, char O[4]);
	void SmallIndType_to_LitEnd_Char(const SmallIndType& I, char O[4]);
	void SmallIndType_to_BigEnd_Char(const SmallIndType& I, char O[4]);
	void BigEnd_Char8_to_PointType(const bool&, char I[8], PointType& O);
	void BigEnd_Char4_to_SmallIndType(const bool&, char I[4], SmallIndType& O);
	void BigEnd_Char4_to_VtxIndType(const bool&, char I[4], VtxIndType& O);
};

/***************************************************************************************/
/* constructor */
MESH_IO::MESH_IO()
{
}

/***************************************************************************************/
/* DE-structor */
MESH_IO::~MESH_IO()
{
}

/***************************************************************************************/
/* write simplex mesh (cells and vertices) to VTK file. (ASCII format) */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MESH_IO::Write_ASCII_VTK(std::string output_filename,
                              std::string vtk_title,
				              const Mesh<CELL_DIM, GEO_DIM>& TM)
{
	// get the file writer
	std::ofstream output;
	// open the file
	output.open(output_filename.c_str(), std::ofstream::out); // default is text format
	// error check
	if (!output)
	{
		std::cerr << std::endl;
		std::cerr << "Write_ASCII_VTK - fatal error!" << std::endl;
		std::cerr << "      Could not open the output file:  '" << output_filename << "'." << std::endl;
		std::exit( 1 );
	}

	// write the header
	Write_VTK_Hdr(output, vtk_title, "ASCII");

	// write vertex coordinates
	const VtxIndType NP = TM.Num_Points();
	Write_VTK_Vtx_Hdr(output, NP, "double");
	
	// write the vertex coordinate data (for each vertex)
	const SmallIndType GD = TM.Geo_Dim();
    for (VtxIndType vv = 0; vv < NP; ++vv)
	{
		// write the vertex's coordinates
		for (SmallIndType kk = 0; kk < GD; ++kk)
		{
			const PointType* VC = TM.Get_Point_coord(vv);
			output << std::scientific << std::setw(26) << std::setprecision(17) << VC[kk];
		}
		// any "extra" coordinates are 0.0 (in VTK, points are always in 3-D)
		for (SmallIndType kk = GD; kk < 3; ++kk)
			output << std::scientific << std::setw(26) << std::setprecision(17) << 0.0;
		// line break
		output << std::endl;
	}

	// blank line
	output << std::endl;
	
	// write the cell data header
	const SmallIndType TD = TM.Top_Dim();
	const CellIndType  NC = TM.Num_Cells();
	const SmallIndType Cell_Order = Write_VTK_Cell_Data_Hdr(output, TD, NC);
	
	// write the cell connectivity data (for each cell)
    for (CellIndType cc = 0; cc < NC; ++cc)
	{
		// write "cell order"
		output << Cell_Order << "  ";
		// write the cell's indices
		const VtxIndType* C_vtx = TM.Get_Cell_vtx(cc);
		for (SmallIndType kk = 0; kk < TD; ++kk)
			output << C_vtx[kk] << "  ";
		// last index and line break
		output << C_vtx[TD] << std::endl;
	}
	
	// blank line
	output << std::endl;
	
	// write cell type information
	const SmallIndType Cell_Label = Write_VTK_Cell_Types_Hdr(output, TD, NC);
    for (CellIndType cc = 0; cc < NC; ++cc)
		output << Cell_Label << std::endl; // write label
	
	//// blank line
	//output << std::endl;

	// close the file
	output.close( );
}

/***************************************************************************************/
/* write simplex mesh (cells and vertices) to VTK file. (BINARY format) */
template <SmallIndType CELL_DIM, SmallIndType GEO_DIM>
void MESH_IO::Write_Binary_VTK(std::string output_filename,
                               std::string vtk_title,
				               const Mesh<CELL_DIM, GEO_DIM>& TM)
{
	// get the file writer
	std::ofstream output;
	// open the file
	output.open(output_filename.c_str(), std::ios::binary); // binary format
	// error check
	if (!output)
	{
		std::cerr << std::endl;
		std::cerr << "Write_Binary_VTK - fatal error!" << std::endl;
		std::cerr << "      Could not open the output file:  '" << output_filename << "'." << std::endl;
		std::exit( 1 );
	}

	// write the header
	Write_VTK_Hdr(output, vtk_title, "BINARY");

	// write vertex coordinates
	const VtxIndType NP = TM.Num_Points();
	Write_VTK_Vtx_Hdr(output, NP, "double");
	
	// write the vertex coordinate data (for each vertex)
	const SmallIndType GD = TM.Geo_Dim();
	const PointType ZERO = 0.0;
	char cc_ZERO[8];
	PointType_to_BigEnd_Char(ZERO, cc_ZERO);
	assert(sizeof(PointType)==8*sizeof(char));
    for (VtxIndType vv = 0; vv < NP; ++vv)
	{
		char coord_char[8];
		// write the vertex's coordinates
		for (SmallIndType kk = 0; kk < GD; ++kk)
		{
			const PointType* VC = TM.Get_Point_coord(vv);
			PointType_to_BigEnd_Char(VC[kk], coord_char);
			output.write(coord_char, sizeof(PointType) );
		}
		// any "extra" coordinates are 0.0 (in VTK, points are always in 3-D)
		for (SmallIndType kk = GD; kk < 3; ++kk)
		{
			output.write(cc_ZERO, sizeof(PointType) );
		}
	}

	// blank line
	output << std::endl << std::endl;
	
	// write the cell data header
	const SmallIndType TD = TM.Top_Dim();
	const CellIndType  NC = TM.Num_Cells();
	const SmallIndType Cell_Order = Write_VTK_Cell_Data_Hdr(output, TD, NC);
	
	// write the cell connectivity data (for each cell)
	assert(sizeof(VtxIndType)==4*sizeof(char));
    for (CellIndType cc = 0; cc < NC; ++cc)
	{
		char cell_char[4];
		// write "cell order"
		SmallIndType_to_BigEnd_Char(Cell_Order, cell_char);
		output.write(cell_char, sizeof(SmallIndType) );
		
		// write the cell's indices
		const VtxIndType* C_vtx = TM.Get_Cell_vtx(cc);
		for (SmallIndType kk = 0; kk < TD; ++kk)
		{
			VtxIndType_to_BigEnd_Char(C_vtx[kk], cell_char);
			output.write(cell_char, sizeof(VtxIndType) );
		}
		// last index
		VtxIndType_to_BigEnd_Char(C_vtx[TD], cell_char);
		output.write(cell_char, sizeof(VtxIndType) );
	}
	
	// blank line
	output << std::endl << std::endl;
	
	// write cell type information
	const SmallIndType Cell_Label = Write_VTK_Cell_Types_Hdr(output, TD, NC);
    for (CellIndType cc = 0; cc < NC; ++cc)
	{
		char label_char[4];
		SmallIndType_to_BigEnd_Char(Cell_Label, label_char);
		output.write(label_char, sizeof(SmallIndType) );
	}
	
	// blank line
	output << std::endl;

	// close the file
	output.close( );
}

/***************************************************************************************/
/* scan the cells and vertices of VTK file and store them in vectors. (ASCII format) */
MeshInterface* MESH_IO::Read_ASCII_VTK(const std::string& input_filename)
{
	SmallIndType TopDim = 0;
	SmallIndType GeoDim = 0;
	std::vector<PointType>    Vtx_Coord;
	std::vector<CellIndType>  Cell_Data;
	CellIndType  Num_CL = 0;
	SmallIndType Cell_Order = 0;
	
	// get the file writer
	std::ifstream input;
	// open the file
	input.open(input_filename.c_str(), std::ifstream::in); // default is text format
	// error check
	if (!input)
	{
		std::cerr << std::endl;
		std::cerr << "Read_ASCII_VTK - fatal error!" << std::endl;
		std::cerr << "      Could not open the input file:  '" << input_filename << "'." << std::endl;
		std::exit( 1 );
	}

	std::string line; // useful!
	
	// read through the file until you get thru all the initial info
	const VtxIndType Num_V = Read_VTK_Initial(input, input_filename, false);

	// read all the vertices and determine geometric dimension
	Vtx_Coord.reserve(3*Num_V);
	PointType vc0, vc1, vc2;
	GeoDim = 0; // initialize
    for (VtxIndType vv = 0; vv < Num_V; ++vv)
	{
		// read the next set of vertex coordinates
		std::getline(input, line);
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;

		// parse it!
		std::stringstream parse_line(line.c_str());
		parse_line >> vc0 >> vc1 >> vc2;
		// store them
		Vtx_Coord.push_back(vc0);
		Vtx_Coord.push_back(vc1);
		Vtx_Coord.push_back(vc2);
		
		// determine geometric dimension
		if ( (GeoDim < 1) && (vc0!=0.0) )
			GeoDim = 1;
		
		if ( (GeoDim < 2) && (vc1!=0.0) )
			GeoDim = 2;
		
		if ( (GeoDim < 3) && (vc2!=0.0) )
			GeoDim = 3;
	}
	
	// read cell data header
	Read_VTK_Cell_Data_Hdr(input, input_filename, Cell_Order, Num_CL);
	
	// read all the cell connectivity data
	Cell_Data.reserve(Num_CL * Cell_Order);
	SmallIndType CO = 0;
	VtxIndType vv_ind = 0;
    for (CellIndType cc = 0; cc < Num_CL; ++cc)
	{
		// read the next set of cell data
		std::getline(input, line);
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;
		
		std::stringstream parse_line(line.c_str());
		parse_line >> CO;
		if (CO!=Cell_Order)
		{
			std::cerr << std::endl;
			std::cerr << "Read_ASCII_VTK - fatal error!" << std::endl;
			std::cerr << "      input file:  '" << input_filename <<
						 "' does *not* have the correct Cell_Order for cell data!" << std::endl;
			std::exit( 1 );
		}
		// read cell vertex indices
		for (SmallIndType kk = 0; kk < Cell_Order; ++kk)
		{
			parse_line >> vv_ind;
			Cell_Data.push_back(vv_ind);
		}
	}
	
	// read cell type header
	Read_VTK_Cell_Type_Hdr(input, input_filename, Num_CL);
	
	// read the cell type information to determine topological dimension
	SmallIndType Cell_Label = 0;
	SmallIndType Prev_Label = 0;
    for (CellIndType cc = 0; cc < Num_CL; ++cc)
	{
		// read the next set of cell type info
		std::getline(input, line);
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;
		
		std::stringstream parse_line(line.c_str());
		parse_line >> Cell_Label;
		if ( (Prev_Label!=0) && (Cell_Label!=Prev_Label) )
		{
			std::cerr << std::endl;
			std::cerr << "Read_ASCII_VTK - fatal error!" << std::endl;
			std::cerr << "      input file:  '" << input_filename <<
						 "' does *not* have *all* cells of the same type!" << std::endl;
			std::exit( 1 );
		}
		Prev_Label = Cell_Label;
	}
	
	// decipher Cell_Label
	std::string Cell_Type_str;
	TopDim = Cell_Label_to_Top_Dim(Cell_Label, Cell_Type_str);
	
	// close the file
	input.close( );
	
	// display mesh information
	Display_Mesh_Info(input_filename, Cell_Type_str, TopDim, GeoDim, Num_CL, Num_V);

	return Convert_to_Mesh_Object(TopDim, GeoDim, Num_CL, Num_V, Cell_Data, Vtx_Coord);
}

/***************************************************************************************/
/* scan the cells and vertices of VTK file and store them in vectors. (BINARY format) */
MeshInterface* MESH_IO::Read_BINARY_VTK(const std::string& input_filename)
{
	SmallIndType TopDim = 0;
	SmallIndType GeoDim = 0;
	std::vector<PointType>    Vtx_Coord;
	std::vector<CellIndType>  Cell_Data;
	CellIndType  Num_CL = 0;
	SmallIndType Cell_Order = 0;
	
	// get the file writer
	std::ifstream input;
	// open the file
	input.open(input_filename.c_str(), std::ifstream::in | std::ifstream::binary);
	// error check
	if (!input)
	{
		std::cerr << std::endl;
		std::cerr << "Read_BINARY_VTK - fatal error!" << std::endl;
		std::cerr << "      Could not open the input file:  '" << input_filename << "'." << std::endl;
		std::exit( 1 );
	}
	
	// read through the file until you get thru all the initial info
	const VtxIndType Num_V = Read_VTK_Initial(input, input_filename, true);

	// read all the vertices and determine geometric dimension
	Vtx_Coord.reserve(3*Num_V);
	PointType vc[3];
	GeoDim = 0; // initialize
	const bool Is_LitEnd = System_Is_Little_Endian();
	assert(sizeof(PointType)==8*sizeof(char));
    for (VtxIndType vv = 0; vv < Num_V; ++vv)
	{
		char coord_char[8];
		// read the 3 coordinates (its always 3 in a VTK file)
		for (SmallIndType kk = 0; kk < 3; ++kk)
		{
			input.read(coord_char, sizeof(PointType) );
			BigEnd_Char8_to_PointType(Is_LitEnd, coord_char, vc[kk]);
		}
		
		// store them
		Vtx_Coord.push_back(vc[0]);
		Vtx_Coord.push_back(vc[1]);
		Vtx_Coord.push_back(vc[2]);
		
		// determine geometric dimension
		if ( (GeoDim < 1) && (vc[0]!=0.0) )
			GeoDim = 1;
		
		if ( (GeoDim < 2) && (vc[1]!=0.0) )
			GeoDim = 2;
		
		if ( (GeoDim < 3) && (vc[2]!=0.0) )
			GeoDim = 3;
	}
	
	// read cell data header
	Read_VTK_Cell_Data_Hdr(input, input_filename, Cell_Order, Num_CL);
	
	// read all the cell connectivity data
	Cell_Data.reserve(Num_CL * Cell_Order);
	SmallIndType CO = 0;
	//VtxIndType vv_ind = 0;
	assert(sizeof(VtxIndType)==4*sizeof(char));
    for (CellIndType cc = 0; cc < Num_CL; ++cc)
	{
		char cell_char[4];
		// read "cell order"
		input.read(cell_char, sizeof(SmallIndType) );
		BigEnd_Char4_to_SmallIndType(Is_LitEnd, cell_char, CO);
		if (CO!=Cell_Order)
		{
			std::cerr << std::endl;
			std::cerr << "Read_BINARY_VTK - fatal error!" << std::endl;
			std::cerr << "      input file:  '" << input_filename <<
						 "' does *not* have the correct Cell_Order for cell data!" << std::endl;
			std::exit( 1 );
		}
		
		// read the cell's indices
		VtxIndType C_vtx[4];
		for (SmallIndType kk = 0; kk < Cell_Order-1; ++kk)
		{
			input.read(cell_char, sizeof(VtxIndType) );
			BigEnd_Char4_to_VtxIndType(Is_LitEnd, cell_char, C_vtx[kk]);
		}
		// last index
		input.read(cell_char, sizeof(VtxIndType) );
		BigEnd_Char4_to_VtxIndType(Is_LitEnd, cell_char, C_vtx[Cell_Order-1]);

		// store cell vertex indices
		for (SmallIndType kk = 0; kk < Cell_Order; ++kk)
		{
			Cell_Data.push_back(C_vtx[kk]);
		}
	}
	
	// read cell type header
	Read_VTK_Cell_Type_Hdr(input, input_filename, Num_CL);
	
	// read the cell type information to determine topological dimension
	SmallIndType Cell_Label = 0;
	SmallIndType Prev_Label = 0;
    for (CellIndType cc = 0; cc < Num_CL; ++cc)
	{
		// read the next cell type info
		char cell_char[4];
		input.read(cell_char, sizeof(SmallIndType) );
		BigEnd_Char4_to_SmallIndType(Is_LitEnd, cell_char, Cell_Label);
		
		if ( (Prev_Label!=0) && (Cell_Label!=Prev_Label) )
		{
			std::cerr << std::endl;
			std::cerr << "Read_BINARY_VTK - fatal error!" << std::endl;
			std::cerr << "      input file:  '" << input_filename <<
						 "' does *not* have *all* cells of the same type!" << std::endl;
			std::exit( 1 );
		}
		Prev_Label = Cell_Label;
	}
	
	// decipher Cell_Label
	std::string Cell_Type_str;
	TopDim = Cell_Label_to_Top_Dim(Cell_Label, Cell_Type_str);
	
	// close the file
	input.close( );
	
	// display mesh information
	Display_Mesh_Info(input_filename, Cell_Type_str, TopDim, GeoDim, Num_CL, Num_V);

	return Convert_to_Mesh_Object(TopDim, GeoDim, Num_CL, Num_V, Cell_Data, Vtx_Coord);
}

/***************************************************************************************/
/* convert raw mesh data into Mesh object. */
MeshInterface* MESH_IO::Convert_to_Mesh_Object(const SmallIndType& TopDim,
											   const SmallIndType& GeoDim,
											   const CellIndType&  Num_CL,
											   const VtxIndType&   Num_V,
											   const std::vector<CellIndType>&  Cell_Data,
                                               const std::vector<PointType>&    Vtx_Coord)
{
	// make sure vertex coordinate data is the proper size
	const VtxIndType Vtx_Coord_Size = 3 * Num_V;
	if (Vtx_Coord_Size != (VtxIndType) Vtx_Coord.size())
	{
		std::cerr << std::endl;
		std::cerr << "Convert_to_Mesh_Object - fatal error!" << std::endl;
		std::cerr << "      Raw mesh data is not consistent with the number" << std::endl;
		std::cerr << "      of vertex coordinates found in the VTK file!" << std::endl;
		std::exit( 1 );
	}

	// make sure cell data is the proper size
	const CellIndType Cell_Data_Size = (TopDim+1) * Num_CL;
	if (Cell_Data_Size != (CellIndType) Cell_Data.size())
	{
		std::cerr << std::endl;
		std::cerr << "Convert_to_Mesh_Object - fatal error!" << std::endl;
		std::cerr << "      Raw mesh data is not consistent with the number" << std::endl;
		std::cerr << "      of cells found in the VTK file!" << std::endl;
		std::exit( 1 );
	}
	
	MeshInterface* M_ptr;
	if ( (TopDim==1) && (GeoDim==1) )
		M_ptr = new Mesh<1,1>;
	else if ( (TopDim==1) && (GeoDim==2) )
		M_ptr = new Mesh<1,2>;
	else if ( (TopDim==1) && (GeoDim==3) )
		M_ptr = new Mesh<1,3>;
	else if ( (TopDim==2) && (GeoDim==2) )
		M_ptr = new Mesh<2,2>;
	else if ( (TopDim==2) && (GeoDim==3) )
		M_ptr = new Mesh<2,3>;
	else if ( (TopDim==3) && (GeoDim==3) )
		M_ptr = new Mesh<3,3>;
	else
	{
		std::cerr << std::endl;
		std::cerr << "Convert_to_Mesh_Object - fatal error!" << std::endl;
		std::cerr << "      Mesh Raw Data:  " <<
					 " has invalid topological and/or geometric dimension!" << std::endl;
		std::exit( 1 );	
	}
	
	// fill in the cell data
    M_ptr->Reserve_Cells(Num_CL);
	for (CellIndType ci = 0; ci < Num_CL; ++ci)
    {
		// read cell vertex indices
		VtxIndType c_vtx[4] = {0, 0, 0, 0};
		for (SmallIndType kk = 0; kk < (TopDim+1); ++kk)
			c_vtx[kk] = Cell_Data[(TopDim+1)*ci + kk];
		M_ptr->Append_Cell(c_vtx);
	}
	
	// fill in the vertex coordinate data
	M_ptr->Init_Points(Num_V);
	for (VtxIndType vi = 0; vi < Num_V; ++vi)
    {
		// read vertex coordinates
		PointType coord[3] = {0.0, 0.0, 0.0};
		for (SmallIndType kk = 0; kk < 3; ++kk)
			coord[kk] = Vtx_Coord[3*vi + kk];
		M_ptr->Set_Coord(vi, coord);
	}

	return M_ptr;
}

/***************************************************************************************/
/* simple output */
void MESH_IO::Display_Mesh_Info(const std::string& input_filename, const std::string& Cell_Type_str,
                                const SmallIndType& TopDim, const SmallIndType& GeoDim,
								const CellIndType&  Num_CL, const VtxIndType& Num_V)
{
	std::cout << "Finished reading mesh data from file: '" << input_filename << "'..." << std::endl << std::endl;
	std::cout << "Mesh Information:" << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	std::cout << "Topological Dimension: " << TopDim << std::endl;
	std::cout << "            Cell Type: " << Cell_Type_str << std::endl;
	std::cout << "      Number of Cells: " << Num_CL << std::endl;
	std::cout << "  Geometric Dimension: " << GeoDim << std::endl;
	std::cout << "   Number of Vertices: " << Num_V << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
}

/***************************************************************************************/
/* simple header */
void MESH_IO::Write_VTK_Hdr(std::ofstream& output,
                            const std::string& vtk_title,
							const std::string& vtk_format)
{
	// write the header
	output << "# vtk DataFile Version 3.0" << std::endl;
	output << vtk_title << std::endl;
	output << vtk_format << std::endl;
	output << "DATASET UNSTRUCTURED_GRID" << std::endl;	
}

/***************************************************************************************/
/* simple header */
void MESH_IO::Write_VTK_Vtx_Hdr(std::ofstream& output, const VtxIndType& Num_Pts, const std::string& format)
{
	// write vertex coordinates header
	output << "POINTS  " << Num_Pts << "  " << format << std::endl;
}

/***************************************************************************************/
/* simple header */
SmallIndType MESH_IO::Write_VTK_Cell_Data_Hdr(
                   std::ofstream& output, const SmallIndType& TopDim, const CellIndType& Num_CL)
{
	// write the cell data header
	const SmallIndType Cell_Order = TopDim + 1;
	const CellIndType Cell_Size = Num_CL * ( Cell_Order + 1 ); // VTK needs this
	output << "CELLS  " << Num_CL << "  " << Cell_Size << std::endl;
	
	return Cell_Order;
}

/***************************************************************************************/
/* simple header */
SmallIndType MESH_IO::Write_VTK_Cell_Types_Hdr(
                   std::ofstream& output, const SmallIndType& TopDim, const CellIndType& Num_CL)
{
	// write cell types information
	output << "CELL_TYPES  " << Num_CL << std::endl;

	// Use topological dimension to determine cell label (for a simplex)
	if (TopDim==1)
		return (SmallIndType)  3; // line (segment)
	else if (TopDim==2)
		return (SmallIndType)  5; // triangle
	else if (TopDim==3)
		return (SmallIndType) 10; // tetrahedron
	else
	{
		std::cerr << std::endl;
		std::cerr << "Write_VTK_Cell_Types_Hdr - fatal error!" << std::endl;
		std::cerr << "      Topological Dimension is invalid:  " << TopDim << "!" << std::endl;
		std::exit( 1 );
	}
	
	return NULL_Small;
}

/***************************************************************************************/
/* simple header */
VtxIndType MESH_IO::Read_VTK_Initial(std::ifstream& input, const std::string& input_filename,
                                     const bool& Is_Binary)
{
	VtxIndType Num_V = 0;
	
	// read through the file until you get thru all the initial info
	std::string line;
	while (!input.eof())
	{
		// read the next line
		std::getline(input, line);
		
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;

		if (Is_Binary)
		{
			// verify
			if ( std::strncmp(line.c_str(),"BINARY",6)==0 )
				continue;

			if ( std::strncmp(line.c_str(),"ASCII",5)==0 )
			{
				std::cerr << std::endl;
				std::cerr << "Read_VTK_Initial - fatal error!" << std::endl;
				std::cerr << "      input file:  '" << input_filename << "' is *not* BINARY format!" << std::endl;
				std::cerr << "      And you indicated to use 'binary'." << std::endl;
				std::exit( 1 );
			}
		}
		else
		{
			// verify
			if ( std::strncmp(line.c_str(),"ASCII",5)==0 )
				continue;

			if ( std::strncmp(line.c_str(),"BINARY",6)==0 )
			{
				std::cerr << std::endl;
				std::cerr << "Read_VTK_Initial - fatal error!" << std::endl;
				std::cerr << "      input file:  '" << input_filename << "' is *not* ASCII format!" << std::endl;
				std::cerr << "      And you indicated to NOT use 'binary'." << std::endl;
				std::exit( 1 );
			}
		}
		
		// verify
		if ( std::strncmp(line.c_str(),"DATASET",7)==0 )
		{
			// verify
			std::size_t found = line.find("UNSTRUCTURED_GRID");
			if (found==std::string::npos)
			{
				std::cerr << std::endl;
				std::cerr << "Read_VTK_Initial - fatal error!" << std::endl;
				std::cerr << "      input file:  '" << input_filename << "' is *not* UNSTRUCTURED_GRID format!" << std::endl;
				std::exit( 1 );
			}
			continue;
		}

		// verify
		if ( std::strncmp(line.c_str(),"POINTS",6)==0 )
		{
			// verify
			std::stringstream parse_line(line.c_str());
			std::string POINTS_str;
		    parse_line >> POINTS_str;
			// get the number of vertices
			parse_line >> Num_V;
			std::string float_str;
		    parse_line >> float_str;
			
			// verify
			if (float_str!="double")
			{
				std::cerr << std::endl;
				std::cerr << "Read_VTK_Initial - fatal error!" << std::endl;
				std::cerr << "      input file:  '" << input_filename <<
				             "' does *not* have 'double' as vertex coordinate type!" << std::endl;
				std::exit( 1 );
			}
			break;
		}
	}
	
	return Num_V;
}

/***************************************************************************************/
/* simple header */
void MESH_IO::Read_VTK_Cell_Data_Hdr(std::ifstream& input, const std::string& input_filename,
                                     SmallIndType& Cell_Order, CellIndType& Num_CL)
{
	std::string line;
	while (!input.eof())
	{
		// read the next line
		std::getline(input, line);
		
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;

		// verify
		std::stringstream parse_line(line.c_str());
		std::string CELLS_str;
		parse_line >> CELLS_str;
		if ( std::strncmp(CELLS_str.c_str(),"CELLS",5)==0 )
		{
			// get the number of cells
			CellIndType All_Cell_Size = 0;
			parse_line >> Num_CL >> All_Cell_Size;
			Cell_Order = ((SmallIndType) All_Cell_Size / Num_CL) - 1;
			break;
		}
	}
	
	if (input.eof())
	{
		std::cerr << std::endl;
		std::cerr << "Read_VTK_Cell_Data_Hdr - fatal error!" << std::endl;
		std::cerr << "      input file:  '" << input_filename <<
					 "' does *not* have any CELL data!" << std::endl;
		std::exit( 1 );
	}
}

/***************************************************************************************/
/* simple header */
void MESH_IO::Read_VTK_Cell_Type_Hdr(std::ifstream& input, const std::string& input_filename,
                                     const CellIndType& Num_CL)
{
	std::string line;
	while (!input.eof())
	{
		// read the next line
		std::getline(input, line);
		
		// if it's a comment or empty, then ignore
		if ( (line[0]=='#') || (line.length()==0) )
			continue;

		// verify
		std::stringstream parse_line(line.c_str());
		std::string CELL_TYPES_str;
		parse_line >> CELL_TYPES_str;
		if ( std::strncmp(CELL_TYPES_str.c_str(),"CELL_TYPES",10)==0 )
		{
			// get the number of cells (as a double check)
			SmallIndType Num_CL_CHK = 0;
			parse_line >> Num_CL_CHK;
			if (Num_CL!=Num_CL_CHK)
			{
				std::cerr << std::endl;
				std::cerr << "Read_VTK_Cell_Type_Hdr - fatal error!" << std::endl;
				std::cerr << "      input file:  '" << input_filename <<
							 "' does *not* have consistent number of cells in the CELL_TYPES section!" << std::endl;
				std::exit( 1 );
			}
			break;
		}
	}
}

/***************************************************************************************/
/* convert cell label to topological dimension and cell type. */
SmallIndType MESH_IO::Cell_Label_to_Top_Dim(const SmallIndType& Cell_Label, std::string& Cell_Type) const
{
	// Use cell label to determine topological dimension (for a simplex)
	if (Cell_Label==3)
	{
		Cell_Type = "line-segment";
		return (SmallIndType)  1; // line (segment)
	}
	else if (Cell_Label==5)
	{
		Cell_Type = "triangle";
		return (SmallIndType)  2; // triangle
	}
	else if (Cell_Label==10)
	{
		Cell_Type = "tetrahedron";
		return (SmallIndType)  3; // tetrahedron
	}
	else
	{
		std::cerr << std::endl;
		std::cerr << "Cell_Label_to_Top_Dim - fatal error!" << std::endl;
		std::cerr << "      Cell_Label is invalid:  " << Cell_Label << "!" << std::endl;
		std::exit( 1 );
	}
	
	return NULL_Small;
}

/***************************************************************************************/
/* determine endian-ness. */
bool MESH_IO::System_Is_Little_Endian()
{
	int n = 1;
	// little endian if true
	if(*(char *)&n == 1)
		return true;
	else
		return false;
}

/***************************************************************************************/
/* convert to little-endian in "char" format. */
void MESH_IO::PointTypeFl_to_LitEnd_Char(const PointTypeFl& I0, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(PointTypeFl) );
	// convert float to UINT32
	const unsigned int INPUT = *(unsigned int*)&I0;
	OUTPUT[0] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[3] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert to BIG-endian in "char" format. */
void MESH_IO::PointTypeFl_to_BigEnd_Char(const PointTypeFl& I0, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(PointTypeFl) );
	// convert float to UINT32
	const unsigned int INPUT = *(unsigned int*)&I0;
	OUTPUT[3] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[0] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert to little-endian in "char" format. */
void MESH_IO::PointType_to_LitEnd_Char(const PointType& I0, char OUTPUT[8])
{
	assert( 8*sizeof(char) == sizeof(PointType) );
	// convert double to UINT64
	const long long unsigned int INPUT = *(long long unsigned int*)&I0;
	OUTPUT[0] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[3] = static_cast<char>((INPUT >> 24) & 0xFF);
	OUTPUT[4] = static_cast<char>((INPUT >> 32) & 0xFF);
	OUTPUT[5] = static_cast<char>((INPUT >> 40) & 0xFF);
	OUTPUT[6] = static_cast<char>((INPUT >> 48) & 0xFF);
	OUTPUT[7] = static_cast<char>((INPUT >> 56) & 0xFF);
}

/***************************************************************************************/
/* convert to BIG-endian in "char" format. */
void MESH_IO::PointType_to_BigEnd_Char(const PointType& I0, char OUTPUT[8])
{
	assert( 8*sizeof(char) == sizeof(PointType) );
	// convert double to UINT64
	const long long unsigned int INPUT = *(long long unsigned int*)&I0;
	OUTPUT[7] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[6] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[5] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[4] = static_cast<char>((INPUT >> 24) & 0xFF);
	OUTPUT[3] = static_cast<char>((INPUT >> 32) & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 40) & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 48) & 0xFF);
	OUTPUT[0] = static_cast<char>((INPUT >> 56) & 0xFF);
}

/***************************************************************************************/
/* convert to BIG-endian to LIT-endian format. */
void MESH_IO::BigEnd_Char_to_LitEnd_Char_32bit(const char INPUT[4], char OUTPUT[4])
{
	OUTPUT[3] = INPUT[0];
	OUTPUT[2] = INPUT[1];
	OUTPUT[1] = INPUT[2];
	OUTPUT[0] = INPUT[3];
}

/***************************************************************************************/
/* convert to BIG-endian to LIT-endian format. */
void MESH_IO::BigEnd_Char_to_LitEnd_Char_64bit(const char INPUT[8], char OUTPUT[8])
{
	OUTPUT[7] = INPUT[0];
	OUTPUT[6] = INPUT[1];
	OUTPUT[5] = INPUT[2];
	OUTPUT[4] = INPUT[3];
	OUTPUT[3] = INPUT[4];
	OUTPUT[2] = INPUT[5];
	OUTPUT[1] = INPUT[6];
	OUTPUT[0] = INPUT[7];
}

/***************************************************************************************/
/* convert to little-endian in "char" format. */
void MESH_IO::VtxIndType_to_LitEnd_Char(const VtxIndType& INPUT, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(VtxIndType) );
	OUTPUT[0] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[3] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert to BIG-endian in "char" format. */
void MESH_IO::VtxIndType_to_BigEnd_Char(const VtxIndType& INPUT, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(VtxIndType) );
	OUTPUT[3] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[0] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert to little-endian in "char" format. */
void MESH_IO::SmallIndType_to_LitEnd_Char(const SmallIndType& INPUT, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(SmallIndType) );
	OUTPUT[0] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[3] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert to BIG-endian in "char" format. */
void MESH_IO::SmallIndType_to_BigEnd_Char(const SmallIndType& INPUT, char OUTPUT[4])
{
	assert( 4*sizeof(char) == sizeof(SmallIndType) );
	OUTPUT[3] = static_cast<char>(INPUT & 0xFF);
	OUTPUT[2] = static_cast<char>((INPUT >> 8) & 0xFF);
	OUTPUT[1] = static_cast<char>((INPUT >> 16) & 0xFF);
	OUTPUT[0] = static_cast<char>((INPUT >> 24) & 0xFF);
}

/***************************************************************************************/
/* convert char[8] to PointType, assuming BIG-endian byte order of char[8]. */
void MESH_IO::BigEnd_Char8_to_PointType(const bool& Is_LitEnd, char INPUT[8], PointType& OUTPUT)
{
	assert( 8*sizeof(char) == sizeof(PointType) );
	if (Is_LitEnd) // is this machine little Endian?
	{
		char char_lit[8];
		BigEnd_Char_to_LitEnd_Char_64bit(INPUT, char_lit);
		// convert to PointType
		OUTPUT = *reinterpret_cast<PointType*>(char_lit);
	}
	else // just convert directly to PointType
		OUTPUT = *reinterpret_cast<PointType*>(INPUT);
}

/***************************************************************************************/
/* convert char[4] to SmallIndType, assuming BIG-endian byte order of char[4]. */
void MESH_IO::BigEnd_Char4_to_SmallIndType(const bool& Is_LitEnd, char INPUT[4], SmallIndType& OUTPUT)
{
	if (Is_LitEnd) // is this machine little Endian?
	{
		char char_lit[4];
		BigEnd_Char_to_LitEnd_Char_32bit(INPUT, char_lit);
		// convert to SmallIndType
		OUTPUT = *reinterpret_cast<SmallIndType*>(char_lit);
	}
	else // just convert directly to double
		OUTPUT = *reinterpret_cast<SmallIndType*>(INPUT);
}

/***************************************************************************************/
/* convert char[4] to VtxIndType, assuming BIG-endian byte order of char[4]. */
void MESH_IO::BigEnd_Char4_to_VtxIndType(const bool& Is_LitEnd, char INPUT[4], VtxIndType& OUTPUT)
{
	if (Is_LitEnd) // is this machine little Endian?
	{
		char char_lit[4];
		BigEnd_Char_to_LitEnd_Char_32bit(INPUT, char_lit);
		// convert to SmallIndType
		OUTPUT = *reinterpret_cast<VtxIndType*>(char_lit);
	}
	else // just convert directly to double
		OUTPUT = *reinterpret_cast<VtxIndType*>(INPUT);
}

// SWW: this is all I need for now.  It might be convenient to read/write to other formats...

#undef MESH_IO

/***/
