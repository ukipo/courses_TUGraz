#!/bin/bash

# If there is no argument, print message and exit
if [ $# -lt 1 ]; then
  echo "$0 needs at least one file"
  exit 1
fi

# If argument is not readable file, print message and exit
if [ $i -r $1 ]; then
  echo "File $1 is readable"
else
  echo "File $1 is not readable"
fi

# Get number of lines in file
# If it is less than 2, print message and exit
count=`cat $1 | wc -l`
if [ $count -lt 2 ]; then
  echo "There is less than 2 lines in $1"
  exit 1
else
  echo "There is 2 or more lines in $1"
# Use head command to print first line
  echo "First line:"
  head -1 $1

# Use tail command to print last line
  echo "Last line:"
  tail -1 $1
fi
