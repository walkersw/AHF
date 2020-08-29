#!/bin/sh

# delete all the test executables
for dir in ./*/; do
  echo -----------------------------------------------------------------------
  echo "Go in this sub-dir: " $dir
  echo "   and delete all *.o executables and all *.vtk files."
  cd "$dir"
  rm *.o
  rm *.vtk
  cd ..
done
echo -----------------------------------------------------------------------
echo "All unit test executables deleted."
exit 0
