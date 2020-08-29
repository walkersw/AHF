#!/bin/sh

# delete all the test executables
for dir in ./*/; do
  echo -----------------------------------------------------------------------
  echo "Go in this sub-dir: " $dir
  echo "   and delete all *.o executables."
  cd "$dir"
  rm *.o
  cd ..
done
echo -----------------------------------------------------------------------
echo "All unit test executables deleted."
exit 0
