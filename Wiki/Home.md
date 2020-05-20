Welcome to the AHF wiki!
========================

AHF: A C++ Mesh Class For The Array-based Half-Facet Data Structure
-------------------------------------------------------------------

NOTE: users should email questions or comments to:  walker@math.lsu.edu.

# Introduction

Unstructured meshes provide the scaffold for many kinds of numerical methods, especially Finite Element Methods (FEM).  This software provides a C++ class for managing and manipulating unstructured simplex meshes in dimensions up to 3.  it is based on the Array-based Half-Facet (AHF) data structure, presented here:

```
Dyedov, V.; Ray, N.; Einstein, D.; Jiao, X.; Tautges, T.
"AHF: array-based half-facet data structure for mixed-dimensional and non-manifold meshes,"
Engineering with Computers, Springer London, 2015, 31, 389-404
```

(See also related references.)

# Features

* The mesh class is split into the topological connectivity (the main part) and a simple class for the vertex coordinates (i.e. the geometry).
* The topological part is implemented for simplex meshes of *any* dimension.
* Initial support for multi-dimensional embedded meshes.
