#!/bin/bash

# state which demo executables to erase
Demo_Dirs=("BaseMesh_0D" "BaseMesh_1D" "BaseMesh_2D" "BaseMesh_3D" "Simple_Mesh_2D")

# run all the demos
for dd in {0..4..1}; do
  demo_dir="${Demo_Dirs[$dd]}"
  echo -----------------------------------------------------------------------
  echo "Delete executable in this sub-dir:" $demo_dir
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cd $demo_dir
  make clean
  rc=$?
  cd ..
  if [ $rc -eq 0 ]
  then
    echo "*Successfully* deleted this demo."
  else
    echo "The error code is:"
	echo $rc
    echo "This demo did not get deleted: " $base_filename
	echo -----------------------------------------------------------------------
	exit 1
  fi
done
echo -----------------------------------------------------------------------
echo "All demo executables erased successfully."
exit 0
