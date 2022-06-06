#!/usr/local/bin/python3

#open file        
yeast_chip = open ('/export/home/up0289/courses/scripting_for_biotechnologists/lab07/yeast_chip.tsv', 'r')

#first line is header
firstline = yeast_chip.readline()     
#append length to the first line
firstline_length = firstline + 'length'
#empty list, in which we will write the contentes
content = []

#iterate over file line by line
for line in yeast_chip:
    #split the line into fields
    fields = line.split()
    #length is the 5th field - 4th field
    length_reg = int(fields[4]) - int(fields[3])
    #append the length at the end of the fileds list
    fields.append(length_reg)

    #append the fields of the line into the content
    content.append(fields)

#sort list by the length (12th element)
yeast_chip_sort=sorted(content, key=lambda x: x[11])
print(firstline_length, *yeast_chip_sort, sep="\n")

yeast_chip.close
