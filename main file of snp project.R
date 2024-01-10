#Main file .. ALL CODES FOR WORKING PROJECT put in here ( Both ALLELE FREQ AND VARIANT CLASSES) 


#Plotting HISTOGRAM of Minor Allele Frequency(MAF) of SNPs for each gene  
#Full code from start to end 
#UCSC Browser 
#starting from here

setwd("C:/Users/sany/Dropbox/My PC (LAPTOP-H8AM3DSJ)/Documents/assignment3/UCSC_BROWSER/MAF_HISTOGRAM_MEAN_MEDIAN")
ab<- read.delim2(file="CCR6_UCSC_REAL.txt", header= T, sep = "\t")       #read the txt file in r by using read.delim2 and separated them into table format by using "/t"
View(ab)


# Load the required library (if not already loaded)
# install.packages("dplyr")  # Run this line if you haven't installed the 'dplyr' package
library(dplyr)

# Remove blank rows in a specific column
cleaned_dataframe <- ab %>%
  filter(!is.na(ab$alleleFreqs) & ab$alleleFreqs != "")    #remove blank rows from allelefreq column 

ROW_COUNT1<- nrow(cleaned_dataframe)  # count the rows 
ROW_COUNT1



#SELECT ONLY THE 1000GENOMES ROWS FROM CLEANED_DATAFRAME
selected_dataframe <- cleaned_dataframe[grep("1000GENOMES", cleaned_dataframe$submitters, ignore.case = TRUE), ]   

ROW_COUNT2<- nrow(selected_dataframe)
ROW_COUNT2

#BY USING FOR LOOP, DELETED THE refUCSC data from the alleles and allelefreq and put the remaining values into another column named as "newallelefreq" 
for (i in 1:851) {
  selected_dataframe$newallelefreq[i] = as.numeric(unlist(strsplit(as.character(selected_dataframe$alleleFreqs[i]), ",")))[selected_dataframe$refUCSC[i]!= unlist(strsplit(as.character(selected_dataframe$alleles[i]),","))]
}

View(selected_dataframe)

#remove the 0 and 1 value from the newallelefreq column
data_clean<- selected_dataframe[selected_dataframe$newallelefreq != 0, ]
data_clean2<- data_clean[data_clean$newallelefreq != 1, ]   # this is optional.. 
data_clean3<- data_clean %>%
  filter(data_clean$class %in% c("single"))                            #remove the "single" character from the class column 


ROW_COUNT3<- nrow(data_clean3)                                          #count the rows of the dataframe
ROW_COUNT3


View(data_clean3)


#Creating a new Dataframe which containing values ( 0.0 to 0.01) for RARE SNPs
df_subset <- subset(data_clean3, data_clean3$newallelefreq >= 0.0 & data_clean3$newallelefreq <= 0.01) #selecting specific range of values from 0.0 to 0.01 (this is for RARE SNPs)



hist(df_subset$newallelefreq,
     labels = paste0(round(hist(df_subset$newallelefreq, plot = FALSE)$counts / length(df_subset$newallelefreq) * 100, 1), "%"),      
     xlim = c(0,0.01), ylim = c(0,200),breaks= 10, xlab = "Frequency", ylab = "No. of SNPs", main = "Histogram of rare SNPs")              #create histogram of newallelefreq data and plot them where the X axis named as freq and the Y axis named as no. of SNPs .. and every single label has the percentage value 
mean_value= mean(df_subset$newallelefreq) #mean
mean_value
median_value= median(df_subset$newallelefreq)  #median
median_value

abline(v = mean_value, col = "red", lwd = 2)  # Creating Vertical line on the plot for mean
abline(v = median_value, col = "blue",  lwd = 2)  # Creating Vertical line on the plot for median
legend("topright", legend = c("Mean", "Median"),
       col = c("red", "blue"), lwd = 2)        #this code is for mentioning which color indicates mean and which indicates median 


#Creating a new dataframe containing values (0.01 to 1.0) for common SNPs

df_subset2 <- subset(data_clean3, data_clean3$newallelefreq >= 0.01 & data_clean3$newallelefreq <= 1.0) #selecting specific range of values from 0.01 to 1.0 (This is for COMMON SNPs)
hist(df_subset2$newallelefreq,
     labels = paste0(round(hist(df_subset2$newallelefreq, plot = FALSE)$counts / length(df_subset2$newallelefreq) * 100, 1), "%"),
     xlim = c(0.01, 1.0), ylim = c(0,20),breaks = 5, xlab = "Frequency", ylab = "No. of SNPs", main = "Histogram of common SNPs")
mean_value= mean(df_subset2$newallelefreq) #mean
mean_value
median_value= median(df_subset2$newallelefreq)  #median
median_value

abline(v = mean_value, col = "red", lwd = 2)  # Vertical line for mean
abline(v = median_value, col = "blue",  lwd = 2)  # Vertical line for median
legend("topright", legend = c("Mean", "Median"),
       col = c("red", "blue"), lwd = 2)






#Plotting barplot of variant classes of SNPs for each genes 

#UCSC genome browser er kaj 
getwd()
setwd("C:/Users/sany/Dropbox/My PC (LAPTOP-H8AM3DSJ)/Documents/assignment3/UCSC_BROWSER/MAF_HISTOGRAM_MEAN_MEDIAN")      #set the working directory where the working file is present 

abcdE<- read.delim2(file= "CXCR3_UCSC_REAL.txt", header= T, sep = "\t")     # read the file by using read.delim2 function
View(abcdE)          # view the file as table format 
abcdE1<- strsplit(abcdE$func, ",")  #split the func column as ","  by using strsplit function and put the value into abcdE1 vector
View(abcdE1) #view the vector
ABCDE2<-table(unlist(abcdE1)) #after splitting , unlist or organized the func column 
par(mar=c(7,4,3,10)) # set the parameter by the help of marzine

barplot(ABCDE2,xlim =c(0,9),ylim = c(0,300), ylab = "No. of SNPs", main = "Barplot of Variant Classes", las=2,offset =0,cex.axis = par("cex.axis"), cex.names = par("cex.axis"))  # create barplot of that unlist data , where all the variant classes are present
