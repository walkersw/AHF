CC=g++
# INC_DIR = ../Demo
CFLAGS=-Wall -Wextra -Werror -g -I.

all: BasePtCoord_2D_example.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

clean:
	rm -rf *o all

# g++ BasePtCoord_2D_example.cpp -o BasePtCoord_2D_example
# g++ -Wall -Wextra -Werror -g BasePtCoord_2D_example.cpp -o BasePtCoord_2D_example
