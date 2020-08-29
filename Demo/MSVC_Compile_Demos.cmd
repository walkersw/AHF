rem script to compile all of the demos.
@echo off

@setlocal EnableDelayedExpansion

set  demo_dirs[0]=BaseMesh_0D
set demo_files[0]=BaseMesh_0D_example

set  demo_dirs[1]=BaseMesh_1D
set demo_files[1]=BaseMesh_1D_example

set  demo_dirs[2]=BaseMesh_2D
set demo_files[2]=BaseMesh_2D_example

set  demo_dirs[3]=BaseMesh_3D
set demo_files[3]=BaseMesh_3D_example

set  demo_dirs[4]=BasePtCoord_2D
set demo_files[4]=BasePtCoord_2D_example

set  demo_dirs[5]=Simple_Mesh_2D
set demo_files[5]=mesh_demo

rem echo %demo_dirs[0]%

@echo on

rem loop through all demo directories
@for /l %%n in (0,1,5) do @(
  @echo --------------------------------------------------
  @echo Compile this directory:  !demo_dirs[%%n]!
  @echo Compile this file:  !demo_files[%%n]!
  cd !demo_dirs[%%n]!
  call ../MSVC_Compile !demo_files[%%n]!
  cd ..
)
@echo --------------------------------------------------
