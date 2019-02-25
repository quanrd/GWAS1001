
read_file <- function(file_path){
  if (!is.null(file_path)) {
    if (file.exists(file_path)) {
      if (endsWith(file_path, ".csv")) {
        dat <- try(read.csv(file_path, header = TRUE, stringsAsFactors = FALSE))
      } else if (endsWith(file_path, ".txt")) {
        dat <- try(read.table(file_path, header = TRUE, stringsAsFactors = FALSE))
      }
      return(dat)
    } else{
      print(paste(file_path, " does not exists.", sep = ""))
    }
  }
  
  return(NULL)
}