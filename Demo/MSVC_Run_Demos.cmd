rem script to run the demos.
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

set  demo_dirs[4]=Simple_Mesh_2D
set demo_files[4]=mesh_demo

rem echo %demo_dirs[0]%

@echo on

rem loop through all demo directories
@for /l %%n in (0,1,4) do @(
  @echo --------------------------------------------------
  @echo Run demo in this directory:  !demo_dirs[%%n]!
  cd !demo_dirs[%%n]!
  call !demo_files[%%n]!.exe
  @if !errorlevel! equ 0 (
    @echo *Successfully* completed this demo: !demo_dirs[%%n]!
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This demo *failed*: !demo_dirs[%%n]!
	@echo --------------------------------------------------
	exit /b
  )
  cd ..
)
@echo --------------------------------------------------
@echo All demos completed successfully.
