
  Implementation of Array Based Half-Facet (AHF) Mesh Data Structure
	                (for simplex meshes only)
                 (c) 08/28/2020, Shawn W. Walker

This code is open source under the BSD license.

DESCRIPTION
========================================================================

This is a C++ implementation of the Array Based Half-Facet (AHF) data structure for unstructured simplicial meshes.  The base class is implemented for *any* dimension through C++ templating.

Source code can be found at:   http://github.com/walkersw/AHF

REMARK: this class uses standard 0-based indexing (C-style).


INSTALLATION AND EXTERNAL LIBRARIES
========================================================================

- Extract the given .zip file.

- You need a C++ compiler.  You can either use MS Visual C++ express edition (in Windows) or the gcc or g++ compiler (in LINUX).  Instructions for both cases are given below.

- AHF now uses the Eigen library found here:  http://eigen.tuxfamily.org/

AHF has been tested with version 3.3.7 of Eigen (probably earlier versions will work).

Therefore, you must *include* the Eigen library when compiling AHF.  Anyone who is willing to use a complex C++ mesh class should be "smart" enough to include an external library.  Nevertheless, I have included some examples below to illustrate.


USAGE
========================================================================

MS Visual C++ (MSVC):
-------------------------

1. See the file "AHF_cmd_example.bat" in the main AHF directory.  It shows how to initialize the MSVC variables and include the AHF source code, as well as the Eigen library.

2. Batch commands are included in the "Unit_Test_src", "Unit_Test_io", and "Demo" sub-directories.

To compile and run the examples, do the following:

Unit Tests: Open a command prompt in the "Unit_Test_src" (or "Unit_Test_io") sub-directory and setup your environment variables and include directories; see the file "AHF_cmd_example.bat" (you will need to modify to suit your setup).

- Next run "MSVC_Compile_Unit_Tests" at the command prompt.
- Then run "MSVC_Run_Unit_Tests" to execute the generated exe's.
- You can run "MSVC_Clean_Unit_Tests" to delete the exe's.
(Note: the unit test files and executables are in sub-dirs of "Unit_Test_src".)

The actual compiling commands can be found in "MSVC_Compile.cmd".

Demos: Open a command prompt in the "Demo" sub-dir and setup your environment variables and include directories as above.

- Next run "MSVC_Compile_Demos" at the command prompt.
- Then run "MSVC_Run_Demos" to execute the generated exe's.
- You can run "MSVC_Clean_Demos" to delete the exe's.
(Note: the demo files and executables are in sub-dirs of "Demo".)


LINUX:
-------------------------

1. See the file "linux_makefile_hdr" in the main AHF directory.  It contains the default compiler options for all the makefile's.  In particular, it contains this line:

EIGEN_DIR=$$HOME/eigen-3.3.7/

which you must modify to fit your installation of the Eigen library.  Nothing else should change.

2. Script files are included in the "Unit_Test_src", "Unit_Test_io", and "Demo" sub-directories.

To compile and run the examples, do the following:

Unit Tests: Open a terminal window in the "Unit_Test_src" (or "Unit_Test_io") sub-directory (your environment variables should already be set up).

- Compile: run "make" at the bash prompt.
- Run: run "./Linux_Run_Unit_Tests.sh" to execute the object files. (make sure you have execute permissions enabled for this script file.)
- You can run "./Linux_Clean_Unit_Tests.sh" to delete the object files.
(Note: the unit test files and object files are in sub-dirs of "Unit_Test_src".)

The actual compiling commands can be found in "makefile" and "linux_makefile_hdr".

Demos: Open a terminal window in the "Demo" sub-directory.

- Compile: run "make" at the bash prompt.
- Run: run "./Linux_Run_Demos.sh" to execute the object files. (make sure you have execute permissions enabled for this script file.)
- You can run "./Linux_Clean_Demos.sh" to delete the object files.
(Note: the demo files and executables are in sub-dirs of "Demo".)


COMPATIBILITY NOTES
========================================================================
This code was developed with MS Visual C++ 2015 and g++, using Notepad++.

Tested on these systems:

-- Windows 10, 64-bit

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

