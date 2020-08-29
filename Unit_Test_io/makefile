include ../linux_makefile_hdr

all: ./Write_Read_Mesh_2D_Nonmanifold/test_WriteRead_2D_Nonmanifold_1.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

# this clean doesn't work
clean:
	rm -rf *o all
    

# g++ test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
# g++ -Wall -Wextra -Werror -g test_Vertex2Halffacet_Class.cpp -o test_Vertex2Halffacet_Class
