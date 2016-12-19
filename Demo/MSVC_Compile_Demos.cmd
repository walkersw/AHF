rem script to compile all of the demos.

rem @echo off

@setlocal EnableDelayedExpansion

rem loop through all demo directories
@for %%i in (Simple_Mesh_2D) do @(
  @echo --------------------------------------------------
  cd %%i
  call ../MSVC_Compile mesh_demo
  cd ..
)
@echo --------------------------------------------------
