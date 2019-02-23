
rm(list = ls())

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 2 & file.exists(args[1])) {
  dat <- try(read.csv(file = args[1], header = TRUE, stringsAsFactors = FALSE))
  
  dat[,1] <- paste("A", dat[,1], sep = "")
  
  write.csv(x = dat, file = args[2], row.names = FALSE, na = "")
} else{
  print("Invalid parameters!!!")
}
