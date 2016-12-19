CC=g++
# INC_DIR = ../Demo
CFLAGS=-Wall -Wextra -Werror -g -I.

all: mesh_demo.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

clean:
	rm -rf *o all

# g++ mesh_demo.cpp -o mesh_demo
# g++ -Wall -Wextra -Werror -g mesh_demo.cpp -o mesh_demo
