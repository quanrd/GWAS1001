
## generate BLUP function
generate_BLUP <- function(dat, by_column = c(1, 2), start_column = 3){
  
  # Convert first column to character
  dat[,1] <- as.character(dat[,1])
  
  # Convert the rest columns to numeric
  for (i in 2:ncol(dat)) {
    dat[,i] <- as.numeric(dat[,i])
  }
  
  # Create lmer formula
  if (length(by_column) > 0) {
    termlabels <- c()
    for (i in 1:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }
  
  # fit the model
  BLUP_out <- apply(dat[,start_column:ncol(dat)], 2, FUN = function(x){
    lme <- lmer(formula = reformulate(termlabels = termlabels, response = "x"), data = dat, REML=TRUE)
    
    # estimate BLUP
    modelblup <- ranef(lme)
    
    # extract BLUP and add grand mean when only one repetition present
    model_line_blup <- modelblup[[1]] + summary(lme)$coefficients[1]
  })
  
  if(length(BLUP_out) > 0){
    for (i in 1:length(BLUP_out)) {
      temp <- as.data.frame(BLUP_out[[i]])
      temp$New <- row.names(temp)
      temp$id <- names(BLUP_out)[i]
      temp <- temp[,c(3, 2, 1)]
      colnames(temp) <- c("id", "Line", "Intercept")
      row.names(temp) <- seq(from = 1, to = nrow(temp), by = 1)
      BLUP_out[[i]] <- temp
      
      if(i==1){
        BLUP_out_df <- BLUP_out[[i]]
      } else{
        BLUP_out_df <- rbind(BLUP_out_df, BLUP_out[[i]])
      }
    }
    
    blup <- dcast(BLUP_out_df, Line ~ id, value.var="Intercept")
    
    colnames(blup)[1] <- colnames(dat)[1]
    blup <- blup[order(as.numeric(gsub("[[:alpha:]]", "", blup[,1]))),]
    row.names(blup) <- seq(from = 1, to = nrow(blup), by = 1)
    
  }
  
  if(exists("blup")){
    return(list("BLUP" = blup))
  } else{
    return(-1)
  }
  
}