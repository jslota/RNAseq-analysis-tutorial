###MMID coding workshop
#
#RNA seq data analysis
#
#Script 1
#Raw read count processing
#
#2022-03-09
#Jessy Slota

###Collect all raw read count files and merge into one matrix
data_files <- Sys.glob("raw read count files/*.tabular") #store paths for all raw read count files
tmp <- list() #create an empty list to store each file
for (i in data_files) { #for loop to load each individual read count file
  X <- gsub(".tabular.*", "", gsub(".*raw read count files/", "", i)) #extract sample name from file name and store in "X"
  tmp[[X]] <- read.delim(i, row.names = 1, header = FALSE) #load read count file "i"
  colnames(tmp[[X]]) <- X #rename column with sample name
  print(X) #print sample name to track progress in console
}
read_counts <- do.call(cbind, tmp) #do.call function collapses all objects within list into one data frame

#Clean up read count file
read_counts <- read_counts[rowMeans(read_counts)>0,] # remove all transcripts that were not detected
read_counts <- read_counts[order(rowMeans(read_counts), decreasing = TRUE),] # order transcripts based on average raw read count

#Save files for further analysis
if (dir.exists("raw data") == FALSE) { dir.create("raw data") } 
write.csv(read_counts, "raw data/raw_read_counts.csv")

