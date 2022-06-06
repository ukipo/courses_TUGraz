#!/bin/bash

# If there is no argument, print a message and exit
if [ $# -lt 1 ]; then
  echo "$0 needs at least one file"
  exit 1
else
	echo "We have at least one file"
fi

# Print a message to say how many files will be copied
echo "Total files: $#"

# Start a for loop to iterate over the arguments,
# copying each file to a backup file 
for i in $@
do
  cp $i backups/$i.bak
done

