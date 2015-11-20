rem script to compile all of the unit tests.

rem @echo off

call MSVC_Compile_AHF test_vertex2halffacet_class
call MSVC_Compile_AHF test_2d_manifold_mesh_1
call MSVC_Compile_AHF test_2d_non_manifold_1
call MSVC_Compile_AHF test_2d_non_manifold_2
