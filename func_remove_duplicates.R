
remove_duplicates <- function(dat, by_column = c(1,2)){
  
  if (exists("temp")) {
    rm(temp)
  }
  
  temp <- list()
  
  for (i in 1:length(by_column)) {
    temp[[i]] <- dat[,i]
  }
  
  # Aggregate data, calculate means for duplicates
  dat <- aggregate(x = dat, by = temp, FUN = "mean")
  
  # Remove the aggregation grouping
  dat <- dat[, (by_column*-1)]
  
  # Re-arrange row names
  row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)
  
  return(dat)
}