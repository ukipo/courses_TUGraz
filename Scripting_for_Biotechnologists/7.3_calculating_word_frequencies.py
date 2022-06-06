#!/usr/local/bin/python3

#counting word frequencies

import string
from collections import Counter


with open("/export/home/up0289/courses/scripting_for_biotechnologists/lab07/lab1.txt") as lab1:
    contents = lab1.read() 
    lowercs = contents.lower()   #make all lowercase
    words = lowercs.split()      #split into individual words
    #get rid of punctuation, special characters etc
    mytable = str.maketrans("", "", string.punctuation) #make mapping table for translating
    justwords = [i.translate(mytable) for i in words] #translate every word 
    
sorted_words = sorted(justwords, key=justwords.count, reverse=True) #sort cleaned up words on the base of their frequency, most frequent first
lisprint = []    #make empty list that you will print
for i in sorted_words:
    if sorted_words.count(i) >= 10:  #if the word appears less than 10 times
        if i not in lisprint:   #and it's not yet in the list to be printed
            lisprint.append(i)  #append it to the list
for i in lisprint:           #for element in the list to be printed
    if i in sorted_words:    #that is also in the soreted words
        print(i, ":", sorted_words.count(i))  #print the element from lisprint and word count from the sorted_words
        
     