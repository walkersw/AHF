include ../linux_makefile_hdr

all: ./Vtx2Halffacet/test_Vertex2Halffacet_Class.o \
./Mesh_2D_Manifold/test_2D_Manifold_Mesh_1.o \
./Mesh_2D_Nonmanifold/test_2D_Nonmanifold_1.o \
./Mesh_2D_Nonmanifold/test_2D_Nonmanifold_2.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

# this clean doesn't work
clean:
	rm -rf *o all
    

# g++ test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
# g++ -Wall -Wextra -Werror -g test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
