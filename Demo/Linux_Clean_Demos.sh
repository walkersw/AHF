#!/bin/sh

# delete all the demo executables
for dir in ./*/; do
  echo -----------------------------------------------------------------------
  echo "Go in this sub-dir: " $dir
  echo "   and delete all *.o executables."
  cd "$dir"
  rm *.o
  cd ..
done
echo -----------------------------------------------------------------------
echo "All demo executables deleted."
exit 0
