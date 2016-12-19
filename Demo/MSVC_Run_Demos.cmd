rem script to run the demos.

@echo off

@setlocal EnableDelayedExpansion

rem loop through all demo directories
@for %%i in (Simple_Mesh_2D) do @(
  @echo --------------------------------------------------
  cd %%i
  call mesh_demo.exe 
  @if !errorlevel! equ 0 (
    @echo *Successfully* completed this demo: %%i
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This demo *failed*: %%i
	@echo --------------------------------------------------
	exit /b
  )
  cd ..
)
@echo --------------------------------------------------
@echo All demos completed successfully.
