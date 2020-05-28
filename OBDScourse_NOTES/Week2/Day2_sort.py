# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
list = [5,8,6,9,1,7,3,2,4]
#Selecting the first as the max

#Question is to find the maximum value across the lsit of numbers
#For each position in the list; 
#run a loop that runs for each position in the list
#Find the maximum and the position of the number
#Swap the number in next position 
#Ultimately have a sorted list
#Step 1 - Find max across the values
#Step 2 - Swap it with the value at last pos
#Step 3 - definition the find max definition1

#line 1: says
#this first definition finds the maximum in the list and the position
# sets the max number to the number in 0 position
#Then sets the max_pos to mean the position of the number 
#Then the loop starts: for each number in the list 
#Add one 
#If the number is more than the maximum (previously set) then
#change the max value to number
#Change the max_pos to the position of new number
#Now return using a tuple where you have to use () to get the number that is the max and it's pos
#print tup1 
#Use definition to find the length
def find_max(list):
    
    max = list[0]
    i = 0    
    max_pos = i
    for number in list[1:]:
        i +=1
        if number > max:
            max = number
            max_pos = i
        #return tuple containing the maxiumum number and the position
        tup1 = (max, max_pos)
        return tup1 
print(find_max(list))

#find length of list
position = len(list)-1

 
#Starting here - next example
#Trying to sort the list now so the numbers are in order I think
#Take the number at <position> in <list> and assign it to <max_pos> in <list>
#Now set maximum to positiion in actuallist   
#
for number in list:
    max,max_pos = find_max(list[0:position+1])
    
    list[max_pos] = list[position]
    list[position] = max
    position = position - 1
    
    print(list)

#I dnt know why my code doens't work


#Now Charlotte is going to git add/commit this script

#Charlotte's code from github - in sort.py
list = [5,8,6,9,1,7,3,2,4]

#run a loop that runs for each position in the list
#find the maximum and end position
#have to first define max and i before the loop
def find_max(list):
    max = list[0]
    i = 0
    max_pos = i
    for number in list[1:]:
        i +=1
        if number > max:
            max = number
            max_pos = i
    tup1 = (max, max_pos)
    return tup1

print(find_max(list))


position = len(list)-1

for number in list:
    max,max_pos = find_max(list[0:position+1])
    
    # Take the number at <position> in <list> and assign it to <max_pos> in <list>
    list[max_pos] = list[position]  
    list[position] = max
    position = position - 1
    
    print(list)


#Moving on - Charlotte and Piyush are making a new branch in terminal
#git checkout -b "new branch"
#git add then git commit 


#Now doing bubble sort
list = [5,8,6,9,1,7,3,2,4]

#First want to look at 1st and 2nd numbers
#if 1 is higher than 2 then swap
#move on and compare 2 and 3, if not higher dont swap

#My try
for i in range(len(list)-1):
    i = 0
    list[position] = list[0]
    if list[i] < list[i+1]:
        list[position] = list[position - 1]
        i +=1 
    elif list[i] > list[i+1]:
        i+=1
    print(list)  
        i +=1
print(list)
    
#Dharam's solution
#So for first bit, say for position i in list
#if the number in i is bigger than i+1
#switch it 
#otherwise keep it the same

for i in range(len(list)-1):
    if list[i]> list[i+1]:
        list[i], list[i+1] = list[i+1], list[i]
    elif list[i]< list[i+1]: 
        list[i], list[i+1] = list[i], list[i+1]
    print(list)

#Need the for within a for to create a loop so it repeats
#This has worked nicely
for j in range (len(list)-1):
    for i in range(len(list)-1):
        if list[i]> list[i+1]:
            list[i], list[i+1] = list[i+1], list[i]
        elif list[i]< list[i+1]: 
            list[i], list[i+1] = list[i], list[i+1]
    print(list)
#Actally dont even need the elif thing
list = [5,8,6,9,1,7,3,2,4]
for j in range (len(list)-1):
    for i in range(len(list)-1):
        if list[i]> list[i+1]:
            list[i], list[i+1] = list[i+1], list[i]
        
    print(list)


list = [5,8,6,9,1,7,3,2,4]
sorted([list])


#open test file using with open
#Next problem
#We want to count GC content and creat an algorith
#For every character, we want to count it 
#We want to work out it's a G or a C
#If it's a GC we count it as a second feature
#Then we do GC count/total count

gc = 0
len_seq = 0
with open("/Users/rhodgson/dnasample.txt", "r") as f:
    for line in f:
        for character in line:
            if character == "A" or character == "G" or character == "C" or character == "T":
                len_seq += 1
                if character == "G" or character == "C":
                    gc += 1
                    
                    
#These print statements need to be after the end of the for loops so they print out after it's finished
print(gc & len_seq)
#Calculate the percentage, a // gives you 0 as an integer         
print(gc/len_seq*100)
                
            



