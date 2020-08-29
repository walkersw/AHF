rem script to delete the demo executables.
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
  @echo Delete demo .exe in this directory:  !demo_dirs[%%n]!
  cd !demo_dirs[%%n]!
  call del !demo_files[%%n]!.exe
  @if !errorlevel! equ 0 (
    @echo *Successfully* deleted this demo: !demo_files[%%n]!
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This demo did not get deleted: !demo_files[%%n]!
	@echo --------------------------------------------------
	exit /b
  )
  cd ..
)
@echo --------------------------------------------------
@echo All demo executables erased successfully.
