rem script to compile all of the unit tests.
@echo off

@setlocal EnableDelayedExpansion

set  unit_test_dirs[0]=Write_Read_Mesh_2D_Nonmanifold
set unit_test_files[0]=test_WriteRead_2D_Nonmanifold_1

rem echo %unit_test_dirs[0]%

@echo on

rem loop through all unit test directories
@for /l %%n in (0,1,0) do @(
  @echo --------------------------------------------------
  @echo Compile in this directory:  !unit_test_dirs[%%n]!
  @echo Compile this file:  !unit_test_files[%%n]!
  cd !unit_test_dirs[%%n]!
  call ../MSVC_Compile !unit_test_files[%%n]!
  cd ..
)
@echo --------------------------------------------------
