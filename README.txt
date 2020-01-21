
  Implementation of Array Based Half-Facet (AHF) Mesh Data Structure
	         (for 1-D, 2-D, 3-D, simplex meshes only)
                 (C) 12/20/2016, Shawn W. Walker

This code is open source under the BSD license.

DESCRIPTION
========================================================================

This is a C++ implementation of the Array Based Half-Facet (AHF) data structure for unstructured simplicial meshes.  The base class is implemented for *any* dimension through C++ templating.

Source code can be found at:   http://github.com/walkersw/AHF

REMARK: this class uses standard 0-based indexing (C-style).

USAGE
========================================================================

1. Extract the given .zip file.

2. You need a C++ compiler.  I used MS Visual C++, express edition.  If you use LINUX, then the gcc or g++ compiler should be fine.

3. If you use MSVC C++, batch commands are included in the "Unit_Test_src", "Unit_Test_io", and "Demo" sub-directories.

To compile and run the .exe's, do the following:

3(a): Open a command prompt in a "Unit_Test_xxx" sub-directory of interest and setup your environment variables by: running "MSVC_Init_Vars" at the command prompt. (Note: you will need to modify "MSVC_Init_Vars.cmd" to suit your setup.)

- Next run "MSVC_Compile_Unit_Tests" at the command prompt.
- Then run "MSVC_Run_Unit_Tests" to execute the generated exe's.

NOTE: the actual compiling commands can be found in "MSVC_Compile.cmd".

3(b): Open a command prompt in the "Demo" sub-dir and setup your environment variables by:  running "MSVC_Init_Vars" at the command prompt.

- Next run "MSVC_Compile_Demos" at the command prompt.
- Then run "MSVC_Run_Demos" to execute the generated exe's.
(Note: the demo files and executables are in sub-dirs of "Demo".)

4. If you use g++ in a Linux environment, then some script files are included in the "Unit_Test_src", "Unit_Test_io", and "Demo" sub-directories.

To compile and run the object files, do the following:

4(a): Open a terminal window in a "Unit_Test_xxx" sub-directory of interest (note: your environment variables should already be set up!).

- Compile: run "make" at the bash prompt.
- Run: run "./Linux_Run_Unit_Tests.sh" to execute the object files. (make sure you have execute permissions enabled for this script file.)

The actual compiling commands can be found in "makefile".

4(b): Open a terminal window in the "Demo" sub-directory of interest.

- Compile: run "./Linux_Make_Demos.sh" at the bash prompt. (make sure you have execute permissions enabled for this script file.)
- Run: run "./Linux_Run_Demos.sh" to execute the object files. (make sure you have execute permissions enabled for this script file.)
(Note: the demo files and executables are in sub-dirs of "Demo".)


COMPATIBILITY NOTES
========================================================================
This code was developed with MS Visual C++ 2015 and g++, using Notepad++.

Tested on these systems:

-- Windows 8, 64-bit

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

