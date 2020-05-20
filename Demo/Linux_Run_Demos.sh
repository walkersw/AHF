#!/bin/bash

# state which demos to run
Demo_Dirs=("BaseMesh_0D" "BaseMesh_1D" "BaseMesh_2D" "BaseMesh_3D" "BasePtCoord_2D" "Simple_Mesh_2D")
Demo_Files=("BaseMesh_0D_example" "BaseMesh_1D_example" "BaseMesh_2D_example" "BaseMesh_3D_example" "BasePtCoord_2D_example" "mesh_demo")

# run all the demos
for dd in {0..5..1}; do
  demo_dir="${Demo_Dirs[$dd]}"
  demo_exe="${Demo_Files[dd]}.o"
  echo -----------------------------------------------------------------------
  echo "Begin running demos in this sub-dir:" $demo_dir
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cd $demo_dir
  ./$demo_exe
  rc=$?
  cd ..
  if [ $rc -eq 0 ]
  then
    echo "*Successfully* completed this demo:" $base_filename
  else
    echo "The error code is:"
	echo $rc
    echo "This demo *failed*:" $base_filename
	echo -----------------------------------------------------------------------
	exit 1
  fi
done
echo -----------------------------------------------------------------------
echo "All demos completed successfully."
exit 0
