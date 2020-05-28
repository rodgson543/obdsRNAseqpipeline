#OBDS Week 4: Base R

#Create character vector with names etc
animals <- c('cat', 'dog', 'giraffe', 'lizard', 'mouse')
animals
#Access element 3
animals['giraffe']
#Replace element 4 with new element 
animals[4] <- 'squirrel'
animals[4]

#Generate an ordered factor of length 6, check levels of the factor 
factor <- factor(c("mouse", "rat", "squirrel", "squirrel", "rat", "mouse"), levels=c("mouse", "rat", "squirrel"), ordered=TRUE)
factor

#Replace element 2 with a value not in levels
factor[2] <- "pigeon"
factor

#Make a list containing character vector from 1, your factor from 2
#Plus a numerical vector and a boolean vector
boolean <- c(TRUE,FALSE,TRUE,FALSE,TRUE)
numerical <- c(1,2,3,4,5,6,2,3,4)
list1 <- list(animals, factor, numerical, boolean)
list1

#Acess differnt elements and then name them usign $
list1[[3]] #Access 3rd object
list1[[3]][1] #Acess first element of 3rd object

#Name objects/elements using $
list1$animals <- list1[[1]] #This adds another element to the list
list1
list1[[5]] <- NULL #Remove object
list1
#Rename objects
names(list1) <- c('animal', 'fact', 'num', 'boo')
list1

#Access using $
list1$animal

levels(list1$fact) <-
  
  
#Create a numeric vector of length 10 
  #
numvect <- c(3,4,15,34,2,6,9,77,13,1)
z <- c(1:10) #Makes a vector of numbers 1 to 10
a <- as.numeric(z) #makes it a num vecotr
a
#Write lapply/sapply to square each eleement in vector
lapply(numvect, function(b) b^2)
class(lapply(numvect, function(b) b^2))
class(sapply(numvect, function(b) b^2))

#Generate a list of length 4 containing both numeric and logical vectors. 
c <- list(a, numvect, c(TRUE, FALSE), z)
c

#Write an lapply or sapply statement to calculate the sum ofthe elements in each vector.
lapply(c, function(f) sum(f))
sapply(c, sum) #get a vector rther rhan list out

#Using your list from 2, write an sapply statement to repeat each
#element of each vector three times Assign the output to a new list.

fish <- sapply(c, function(c) rep(c,each=3)) #repeats 3 times
#Could do: sapply(c, rep, each=3)
fish

#Generate matrix with 5 rows containing numbers 2:100 in increments of 2, by row
matrix <- matrix(seq(2, 100, by=2), nrow=5, byrow = TRUE)
matrix

#means of row using apply 
apply(matrix, 1, mean)
#(and the sum of each column)
apply(matrix,2,sum)

#Second matrix, 10 rows, 6 columns, fill with numbers of whatever
matrix2 <- matrix(1:60, nrow=10, ncol=6, byrow=TRUE)
matrix2
#Calculate transponse of matrix
matrix2t <- t(matrix2)
matrix2t
#Join transposd matrix to matrix from 1(join by row) 
matrix_bind <- rbind(matrix, matrix2t)

#Check dimensions
dim(matrix_bind)

#Convert to dataframe 
DF <- as.data.frame(matrix_bind)
DF
class(DF)

#Data frame exercises 
coding_bed <- read.table("/Users/rhodgson/OBDStestdata/Week5/coding_gene_region.bed",  header=FALSE, sep="\t", stringsAsFactors = FALSE)
#Check the dimensions of the data frame and the class of each variable.
coding_bed
dim(coding_bed)
class(coding_bed)
#Check class of each variable using apply - have to use sapply here as matrix only contain one data type - apply will coerce this
sapply(coding_bed, typeof) #Change columns to be numeric
class(coding_bed)
#Convert
#coding_bed[[2]] <- as.numeric(coding_bed[[2]])
#apply(coding_bed, 2, typeof) #Change columns to be numeric

#Add column names
colnames(coding_bed)

cols <- c("Chr", "Start", "Stop", "Name","Score", "Strand") 
colnames(coding_bed) <- cols

#Add a new column containing the lenght of each genomic interval and
#sort this column from largest to smallest using a base R function

#Length will be stop minus start
coding_bed$Length <- coding_bed$Stop - coding_bed$Start

#Sort (need the #means you dont want to subset) #This has sorted by ascendinh
sorted_bed <- coding_bed[order(coding_bed$Length, decreasing=TRUE), ]
sorted_bed
#Adding a minus dirctly after order will have does this step for you

#Extract element at row 30, col 3
sorted_bed[30, 3]

#ExtrExtract the second column by index and by name (using both [] and $)

sorted_bed[[2]] #by index
sorted_bed$Start #by name

#On which chromosome is the largest interval? Output just the
#chromosome value and store in the variable max_chrom
#Subsetting column that contains chr value, gives you chr value
max_chrom <- coding_bed$Chr[coding_bed$Length == max(coding_bed$Length)]
max_chrom
coding_bed$Length == max(coding_bed$Length) #Logical value where it prints out which value contains the maxlength
#above is the condition, then we subset for the col on the condition here 
            

#7. Subset the data frame to contain only regions with a length from
#100,001-200,000 bp - assign to a new variable. Write your subset
#data frame to a tab separated  le (include column names but not rownames).
             
subset_bed <- sorted_bed[sorted_bed$Length >100000 & sorted_bed$Length <= 200000, ]
View(subset_bed)
#Write to tabsep file (col name but not rownmaes)
write.table(subset_bed, "/Users/rhodgson/OBDStestdata/Week5/subset_bed.txt", sep ="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)


#In the original data frame, replace the score value with 100 for genomic intervals
#on chr4 or chr17 that are on the + strand and longer than 200,000 bp. 
sorted_bed$Score[sorted_bed$Chr == "chr4" | sorted_bed$Chr == "chr17" & sorted_bed$Strand == "+" & sorted_bed$Length > 200000] <- 100

#Count the number of regions that have a score of 100.
length(which(sorted_bed$Score ==100))

#Add row to dataframe
row <- data.frame(Chr="chr6", Start= 1000000, Stop=1000200, Name="ENSG00000564732", Score=100, Strand="+", Length=200)
extend_bed <- rbind(sorted_bed, row)
extend_bed
class(row$Chr)
sapply(extend_bed, typeof) #double is same as int
#Remove score variable from datafram
extend_bed$Score <- NULL
tail(extend_bed)
#Can also use extend_bed[,-5] or drop(sorted_bed)
#Use apply to find max of each column
apply(extend_bed, 2, max)
      
#Loops etc
colours_vector <- c("red", "orange", "purple", "yellow", "pink", "blue")
#Write a loop to print the colours in colours_vector with four characters
for (i in colours_vector){
  if (nchar(i) == 4){
    print(i)}
}

#Write a loop to print out the colours at even positions of the
#colours_vector (loop should work for a vector of any length)
for (i in seq_along(colours_vector)){
  if (i %% 2 == 0){
  print(colours_vector[i])}
}
seq_along(colours_vector[i])

#Write a function that uses a for loop to replace all instances of "e" for
#"o" in the colours_vector (note in reality you would do this with sapply)

replace_function <- function(var){
  newcolour <-gsub("e", "o", colours_vector)
  return(newcolour)
}
replace_function()

  
#Write a function that uses a for loop to calculate the mean of a
#numeric vector of any length (use of the mean() function is banned)
vector <- c(1,2,4,5,6,7,8, 436,23,57,25,8547,234,1434)

newmean <- function(var){
  sum = 0
  for (i in var){
    sum = sum+i}
  return(sum/length(var))
}
newmean(vector)

#Advanced: write a function that returns the number of vowels in each
#element of the colours_vector , but only for elements with fewer than
#six characters
nrvowels <-function(myvector){
  vowel = c("a", "e", "o", "i" ,"u")
  output = c()
  for (i in myvector){
    number = 0
    if (nchar(i)< 6){
    split <- unlist(strsplit(i, ""))
    for (x in split){
    if (x %in% vowel){
      number = number+1}}
    output = c(output, number)
    }}
  return(output)}

#Test function
nrvowels(colours_vector)

#Sum 
nrvowels <-function(myvector){
  vowel = c("a", "e", "o", "i" ,"u")
  output = c()
  for (i in myvector){
    number = 0
    if (nchar(i)< 6){
      split <- unlist(strsplit(i, ""))
      for (x in split){
        if (x %in% vowel){
          number = number+1}}
      output = c(output, number)
    }}
  return(sum(output))}
nrvowels(colours_vector)


#Plotting in base R 
data(mtcars)
#Draw a bar plot showing the number of cars with 3, 4 or 5 forward
#gears. Hint: use table() function to get data into the correct format.
#Add x and y axis labels
#Add a main title spread across two lines
#Change the colour and width of the bars
#Add a horizontal line at y = 6
#Change plot margins to 6, 6, 5, 5
x<- table(mtcars$gear) #Could do within function
#par(mar = c(6,6,5,5))
barplot(x, xlab="Gears", ylab = "Frequency", main = "Freq of cars by gears" , col = "green", border = "black", abline(h=6, col="blue"))

#Next
#Generate a scatter plot of mpg vs. hp coloured by gear values
#Change points to  lled circles and increase their size
#Add x and y axis titles and change size
#Change colours to red, green and blue
#Add a legend
plot(mtcars$mpg, mtcars$hp, col = mtcars$gear, xlab = "Miles per gal", ylab="Horsepower", cex )

#Generate a scatter plot of mpg vs. hp coloured by gear values
##Add x and y axis titles (xlab,ylab) and change size (cex.lab)
#Change points to filled circles (pch) and increase their size(cex)
#To change colours to red, green and blue we gsub each value of variable gears with colours. Then use the same in arguments for legend to ensure congruent colour assignments between legend and plot )
gears <- mtcars$gear
gears <- gsub("3", "green", gears)
gears <- gsub("4", "blue", gears)
gears<- gsub("5", "red",gears)
plot(mtcars$mpg,mtcars$hp, col= gears, xlab = "miles per galleon", ylab = "Horsepower",pch=16,cex = 1.5,cex.lab=1.2)
#Add a legend
legend(25, y=350, legend = c("3", "4", "5"), fill = c("blue", "green", "red"))



#Plot the two plots that we have just made side-by-side by changing the global graphical parameters
par()$mfrow #check global parameter mfrow 
par(mfrow = c(1,2)) #change to  2 columns
#Now if we re-plot the 2 graphs, they appear side by side (they appear squished but look fine on zooming)



