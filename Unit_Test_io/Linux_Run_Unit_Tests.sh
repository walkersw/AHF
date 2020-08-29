#!/bin/bash

# state which unit tests to run
Unit_Test_Dirs=("Write_Read_Mesh_2D_Nonmanifold")
Unit_Test_Files=("test_WriteRead_2D_Nonmanifold_1")

# run all the unit tests
for dd in {0..0..1}; do
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
