rem script to run the unit tests.
@echo off

@setlocal EnableDelayedExpansion

set  unit_test_dirs[0]=Write_Read_Mesh_2D_Nonmanifold
set unit_test_files[0]=test_WriteRead_2D_Nonmanifold_1

rem echo %unit_test_dirs[0]%

@echo on

rem loop through all unit test directories
@for /l %%n in (0,1,0) do @(
  @echo --------------------------------------------------
  @echo Run unit test in this directory:  !unit_test_dirs[%%n]!
  cd !unit_test_dirs[%%n]!
  call !unit_test_files[%%n]!.exe
  @if !errorlevel! equ 0 (
    @echo *Successfully* completed this unit test: !unit_test_files[%%n]!
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This unit test *failed*: !unit_test_files[%%n]!
	@echo --------------------------------------------------
	exit /b
  )
  cd ..
)
@echo --------------------------------------------------
@echo All unit tests completed successfully.
