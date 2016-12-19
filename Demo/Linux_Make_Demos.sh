#!/bin/sh

# state which demos to compile
Demos="Simple_Mesh_2D"

# run all the demos
for demo_dir in $Demos; do
  echo -----------------------------------------------------------------------
  echo "Begin compiling demos in this sub-dir:" $demo_dir
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cd $demo_dir
  make
  rc=$?
  cd ..
  if [ $rc -eq 0 ]
  then
    echo "*Successfully* compiled this demo:" $base_filename
  else
    echo "The error code is:"
	echo $rc
    echo "This demo *failed* to compile:" $base_filename
	echo -----------------------------------------------------------------------
	exit 1
  fi
done
echo -----------------------------------------------------------------------
echo "All demos compiled successfully."
exit 0
