#!/bin/bash

# If no input file is given, ask user to give one
if [ $# -lt 1 ]; then
  read -p "Input file? " file
else
  file=$1
fi

# Print header line telling what will be printed:
# Line number, Number of words, Line
echo "line	words	line"

# Initialize line counter variable to 0
k=0



# Read one line of input file, assign it to a variable,
# while there is a line remaining

while read line
do
  # Increase line counter by 1
   ((k++))
   wrds=$(echo "$line" | wc -w)
  # Print line number, number of words and line itself, separate by TAB
   echo "$k	" $wrds "	$line"
done < $file
