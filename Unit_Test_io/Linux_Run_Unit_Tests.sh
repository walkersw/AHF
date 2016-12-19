#!/bin/sh

# state which tests to run
Unit_Tests="./test_write_read_2d_non_manifold_1.o"

# run all the tests
for filename in $Unit_Tests; do
  echo -----------------------------------------------------------------------
  TEMP_NAME=$(basename $filename)
  extension="${TEMP_NAME##*.}"
  base_filename="${TEMP_NAME%.*}"
  echo "Begin running this unit test:" $base_filename
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  $filename
  rc=$?
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
