
  Implementation of Array Based Half-Facet (AHF) Mesh Data Structure
	         (for 1-D, 2-D, 3-D, simplex meshes only)
                 (C) 12/09/2016, Shawn W. Walker

This code is open source under the BSD license.

DESCRIPTION
========================================================================

This is a C++ implementation of the Array Based Half-Facet (AHF) data structure for unstructured simplicial meshes.  The base class is implemented for *any* dimension through C++ templating.

Source code can be found at:   http://github.com/walkersw/AHF

USAGE
========================================================================

1. Extract the given .zip file.

2. You need a C++ compiler.  I used MS Visual C++, express edition.  If you use LINUX, then the gcc or g++ compiler should be fine.

3. If you use MSVC C++, several batch commands are included in the "Unit_Test" sub-directory.

To compile and run the .exe's, do the following:

- Open a command prompt in "./Unit_Test" and setup your environment variables by:

running "MSVC_Init_Vars" at the command prompt. (Note: you will need to modify "MSVC_Init_Vars.cmd" to suit your setup.)

- Next run "MSVC_Compile_AHF_Unit_Tests" at the command prompt.
- Then run "MSVC_Run_AHF_Unit_Tests" to execute the generated exe's.

NOTE: the actual compiling commands can be found in "MSVC_Compile.cmd".

4. If you use g++ in a Linux environment, then some script files are included in the "Unit_Test" sub-directory.

To compile and run the object files, do the following:

- Open a terminal window in "./Unit_Test"; note: your environment variables should already be set up!
- Compile: type "make" at the bash prompt and press enter.
- Run: type "./Linux_Run_AHF_Unit_Tests.sh" to execute the object files. (make sure you have execute permissions enabled for this script.)

The actual compiling commands can be found in "makefile".

5. NOTE: this class uses 0-based indexing!


COMPATIBILITY NOTES
========================================================================
This code was developed with MS Visual C++ 2010 and g++, using Notepad++.

Tested on these systems:

-- Windows 7, 8, 64-bit

-- LINUX KDE/Ubuntu, 64-bit


BUG REPORTS AND FEEDBACK
========================================================================
Please report any problems and/or bugs to:  walker@math.lsu.edu


ACKNOWLEDGEMENTS
========================================================================

This implementation is based on the paper:

Dyedov, V.; Ray, N.; Einstein, D.; Jiao, X.; Tautges, T.
"AHF: array-based half-facet data structure for mixed-dimensional and non-manifold meshes,"
Engineering with Computers, Springer London, 2015, 31, 389-404

