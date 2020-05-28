# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#OBDS course week 2; python introduction

#This will print the numbers between 1 to 9
for i in range(1,10,1):
        print(i)
       
#This will print the square        
for i in range(1,10,1):
        print(i**2)
    
#Now want to print the natural log of 10
import math

math.log(10,2)

math.log1p(1)
math.sqrt(1234)
math.round(0.05)
round(0.05)
math.ceil(0.05)
math.floor(0.05)
#prints 2 to the power of 10
math.pow(2, 10)

#print valuwe of pi
math.pi 
math.e


#print time
import time
help(time)
clock()
import datetime
help(datetime)

now = datetime.datetime.now()

#print date and time wihtout millisecons
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print("date and time =", dt_string)

#now print out just date
dt2 = now.strftime("%d/%m/%Y ")
print("date =", dt2)

#compute number of days since 1st of january
from datetime import date
f_date = date(2020, 1, 1)
l_date = date(2020, 5, 4)
delta = l_date - f_date
print(delta.days)


#Assign"Oxford Biomedical Data science Training programme to a vairbalae"

OBDS = "Oxford Biomedical Data Science Training Programme"
obds= OBDS.lower()
print(obds)
#this splits the string into words
splitOBDS = obds.split()
print(splitOBDS)

#print each word, tab space and it's length
for word in splitOBDS:
    print (word+"\t"+str(len(word)))
    
#print Print the first half of the words in the list.
for word in splitOBDS
    print(f"{word}\t{len(word)}")

#Open file
import os
os.getcwd()
os.chdir('/Users/rhodgson/')

ERR =open('ERR1755082.test.sam', 'r')
#Read everything in the file by line
ERR_data =  ERR.readlines()
print(ERR_data)

#Close the file cause you dont need it as you've already read the data into a new variable 
ERR.close()
#Check how many lines in the file
len(ERR_data)
#print the no.characters. then number of tab-sparated fields
for line in ERR_data:
    print (str(len(line)) + "\t"+ str(line.count("\t")))
# For index starting at 0 and ending at 9
for i in range(0,10,1):
# Get line at that index
for i in range(0,10,1):  
    print(ERR_data[i])
    print (str(len(ERR_data[i])) + "\t"+ str(ERR_data[i].count("\t")))
# print number of character and number of fields
for i in range(0,10,1):  
    #print(ERR_data[i])
    print (str(len(ERR_data[i])) + "\t"+ str(ERR_data[i].count("\t")))
#so the first bit just says for each of the line from 0 to 10
#do the 



ERR2 =open('ERR1755082.test.sam', 'r')
#Read everything in the file by line
ERR2_data =  ERR2.readlines()
print(ERR2_data)

#Close the file cause you dont need it as you've already read the data into a new variable 
ERR2.close()
#Check how many lines in the file
len(ERR2_data)
#print the no.characters. then number of tab-sparated fields
for line in ERR2_data:
    print (str(len(line)) + "\t"+ str(line.count("\t")))
# For index starting at 0 and ending at 9
for i in range(0,10,1):
# Get line at that index
for i in range(0,10,1):  
    print(ERR2_data[i])
    print (str(len(ERR2_data[i])) + "\t"+ str(ERR2_data[i].count("\t")))
# print number of character and number of fields
for i in range(0,10,1):  
    #print(ERR_data[i])
    print (str(len(ERR_data[i])) + "\t"+ str(ERR_data[i].count("\t")))
#so the first bit just says for each of the line from 0 to 10
#do the 

#Import logging 
import logging as L

#Test for the debugger in spyder

items = [1, 2, 3, 'test', 4]
for i in range(len(items)):
    item = items[i]
    value = item // 2

#Now going to have another go at looking at some python problems from coding bat
#Warmup-1 #Sleep_in
#The parameter weekday is True if it is a weekday, and the parameter vacation is True if we are on vacation. 
#We sleep in if it is not a weekday or we're on vacation. Return True if we sleep in.
#def sleep_in(weekday, vacation):
#    if vacation or weekend is true then sleep_in is true:
#        elif weekday is true then sleep_in is false         

#Worked through with youtube tutorial
def sleep_in(weekday, vacation):
    if weekday == False or vacation == True:
        return True
    else:
        return False
#Test this now     
sleep_in(True, False)
sleep_in(False, False)

#Warmup-2 #string_times
#Given a string and a non-negative int n, 
#return a larger string that is n copies of the original string.
#string_times('Hi', 2) → 'HiHi'
#string_times('Hi', 3) → 'HiHiHi'
#string_times('Hi', 1) → 'Hi'



s = "Hi"
print(s+s)

def string_times(str, n):
    return str*n

#Test
string_times('Hi', 4)


#String-1 #Hello_name
#given a string name, return a greeting of the form Hello Bob
#hello_name('Alice') → 'Hello Alice!'

#I did this one myself!!
def hello_name(name):
    return("Hello "+ name)

hello_name("Rose")

#List-1 #first_last6
#Given an array of ints, return True if 6 appears as
#either the first or last element in the array. 
#The array will be length 1 or more.
a= [1,2,3,4,5,6]
a[0]
#print second from last number
len(a)-1
#print last number
a[-1]

#WWoooohooo!
def first_last6(nums):   
    if nums[0] == 6 or nums[-1] == 6:
        return True
    else:
        return False
     
first_last6([1,2,3,4,5,6])


#Next one is from strings_2

#Given a string, return a string where for every char in the original,there are two chars.
#double_char('The') → 'TThhee'
#double_char('AAbb') → 'AAAAbbbb'
#double_char('Hi-There') → 'HHii--TThheerree'
str = "ABC"

def double_char(str):
    for c in str:
        for i in str:
            return i&i
            i+=1
        print(str)
        
        
double_char("Hello")

str = "ABC"

def double_char(str):
    to_return = ''
    for c in str:
        to_return += c*2
    return to_return


double_char("Hello")

#Next problem: List-2 count_evens
#Return the number of even ints in the given array. 
#Note: the % "mod" operator computes the remainder, e.g. 5 % 2 is 1.

#count_evens([2, 1, 2, 3, 4]) → 3
#count_evens([2, 2, 0]) → 3
#count_evens([1, 3, 5]) → 0

def count_evens(nums):
    
    
nums = 
    
    
    
    
    
    
    
    
    
    
    
    
    