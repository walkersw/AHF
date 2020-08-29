include ../linux_makefile_hdr

all: ./BaseMesh_0D/BaseMesh_0D_example.o \
./BaseMesh_1D/BaseMesh_1D_example.o \
./BaseMesh_2D/BaseMesh_2D_example.o \
./BaseMesh_3D/BaseMesh_3D_example.o \
./BasePtCoord_2D/BasePtCoord_2D_example.o \
./Simple_Mesh_2D/mesh_demo.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

# this clean doesn't work
clean:
	rm -rf *o all
    

# g++ test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
# g++ -Wall -Wextra -Werror -g test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
