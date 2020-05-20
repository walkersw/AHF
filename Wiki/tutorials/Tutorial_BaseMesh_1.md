The BaseMesh Class
==================

The core class of this library is the BaseMesh class.  It only contains topological information about the simplex mesh; no geometric coordinate information is stored.  The BaseMesh class is templated on the topological dimension of the mesh, starting at 0.  The dimension can be as large as the largest integer *long int* on your machine!  However, for practical reasons, the dimension should be less than, say, 10.  The library is not optimized when the dimension is very high.

Under the `Demos` directory, go through the demos given in the sub-dirs: `BaseMesh_XD`, where X is the topological dimension.  I suggest the following order:

* BaseMesh_1D
* BaseMesh_0D
* BaseMesh_2D
* BaseMesh_3D

(Starting with 0D is a little too weird.)

The codes are very well commented.  You should also work through the examples by hand to improve your understanding.
