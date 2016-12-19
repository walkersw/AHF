CC=g++
# INC_DIR = ../io_code
CFLAGS=-Wall -Wextra -Werror -g -I.

all: test_write_read_2d_non_manifold_1.o

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
# NOTE: $< refers to the LEFT of :
#       $@ refers to the RIGHT of :

clean:
	rm -rf *o all

# g++ test_write_read_2d_non_manifold_1.cpp -o test_write_read_2d_non_manifold_1
# g++ -Wall -Wextra -Werror -g test_write_read_2d_non_manifold_1.cpp -o test_write_read_2d_non_manifold_1
