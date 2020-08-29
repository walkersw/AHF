#!/bin/bash

# state which unit tests to run
Unit_Test_Dirs=("Vtx2Halffacet" \
"Mesh_2D_Manifold" \
"Mesh_2D_Nonmanifold" \
"Mesh_2D_Nonmanifold")
Unit_Test_Files=("test_Vertex2Halffacet_Class" \
"test_2D_Manifold_Mesh_1" \
"test_2D_Nonmanifold_1" \
"test_2D_Nonmanifold_2")

# run all the unit tests
for dd in {0..3..1}; do
  ut_dir="${Unit_Test_Dirs[$dd]}"
  ut_exe="${Unit_Test_Files[dd]}.o"
  echo -----------------------------------------------------------------------
  echo "Begin running unit tests in this sub-dir:" $ut_dir
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cd $ut_dir
  ./$ut_exe
  rc=$?
  cd ..
  if [ $rc -eq 0 ]
  then
    echo "*Successfully* completed this unit test:" $base_filename
  else
    echo "The error code is:"
	echo $rc
    echo "This unit test *failed*:" $base_filename
	echo -----------------------------------------------------------------------
	exit 1
  fi
done
echo -----------------------------------------------------------------------
echo "All unit tests completed successfully."
exit 0
