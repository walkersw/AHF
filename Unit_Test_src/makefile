CC=g++
# INC_DIR = ../src_code
CFLAGS=-Wall -Wextra -Werror -g -I.

all: test_vertex2halffacet_class.o test_2d_manifold_mesh_1.o test_2d_non_manifold_1.o test_2d_non_manifold_2.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

clean:
	rm -rf *o all

# g++ test_vertex2halffacet_class.cpp -o test_vertex2halffacet_class
# g++ -Wall -Wextra -Werror -g test_vertex2halffacet_class.cpp -o test_vertex2halffacet_class
