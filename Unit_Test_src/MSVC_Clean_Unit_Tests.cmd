rem script to delete the unit test executables.
@echo off

@setlocal EnableDelayedExpansion

set  unit_test_dirs[0]=Vtx2Halffacet
set unit_test_files[0]=test_Vertex2Halffacet_Class

set  unit_test_dirs[1]=Mesh_2D_Manifold
set unit_test_files[1]=test_2D_Manifold_Mesh_1

set  unit_test_dirs[2]=Mesh_2D_Nonmanifold
set unit_test_files[2]=test_2D_Nonmanifold_1

set  unit_test_dirs[3]=Mesh_2D_Nonmanifold
set unit_test_files[3]=test_2D_Nonmanifold_2

rem echo %unit_test_dirs[0]%

@echo on

rem loop through all unit test directories
@for /l %%n in (0,1,3) do @(
  @echo --------------------------------------------------
  @echo Delete unit test .exe in this directory:  !unit_test_dirs[%%n]!
  cd !unit_test_dirs[%%n]!
  call del !unit_test_files[%%n]!.exe
  @if !errorlevel! equ 0 (
    @echo *Successfully* deleted this unit test: !unit_test_files[%%n]!
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This unit test did not get deleted: !unit_test_files[%%n]!
	@echo --------------------------------------------------
	exit /b
  )
  cd ..
)
@echo --------------------------------------------------
@echo All unit test executables erased successfully.
