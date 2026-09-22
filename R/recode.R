#' Encode categorical variables
#'
#' The function is used to encode categorical variables in a data frame, based on the specified encoding type.
#' It supports three encoding types: continuous, dummy, and effect coding.
#' The input data can be in either long or wide format.
#'
#' @name recode
#' @param data a data frame containing the variables to be encoded.
#' The data frame can be in either long or wide format. Variables can be of type numeric or character.
#' @param name a character vector of variable names to encode.
#' @param ctype a character vector specifying the encoding type for each variable:
#' "C" for continuous (no encoding), "D" for dummy coding, and "E" for effect coding.
#' @param reflevel a character vector specifying the reference levels for encoding.
#' For variables of type character, the default `reflevel` is the alphabetically first level.
#' For numeric variables, the default `reflevel` is the level with smallest numeric value.
#' @param ... additional arguments (currently not used).
#'
#' @return A data frame with the encoded variables based on the specified encoding types.
#' @examples
#' data(health, package = "groltir")
#' name <- c("type", "timeSpan", "effectivenss", "adverseEffects")
#' ctype <- c("D", "E", "C", "D")
#' reflevel <- c(1, 1, NA, 2)
#' health_coded <- recode(data = health, name, ctype, reflevel) # data in wide format
#' head(health_coded)
#' @export
#'
recode <- function(data, name, ctype, reflevel, ...){
  if (length(name) != length(ctype)) {
    stop("Error: Length of 'name' and 'ctype' are not consistent.")
  }

  newdf <- data.frame(matrix(nrow = nrow(data), ncol = 0))
  for (i in 1:length(name)){
    var <- name[i]
    search <- paste("^", var, sep='')
    index <- grep(search, names(data))
    rest <- gsub(search, "", names(data[,index]))

    if (ctype[i] == "C"){
      newnames <- colnames(data)[index]
      newdf[,newnames] <- data[,index]
    }else if(ctype[i] == "D"){
      lvls <- sort(unique(unlist(data[,index])))
      if (!is.null(reflevel[i])){
          index_ref <- which(lvls %in% reflevel[i])
          lvls <- c(lvls[-index_ref], lvls[index_ref])}else{
          lvls <- c(lvls[-1], lvls[1])
        }
      for (j in 1:(length(lvls)-1)){
        newnames <- paste(var, "_level", lvls[j], rest, sep='')
        newdf[,newnames] <- sapply(data[,index], function(x) ifelse(x==lvls[j],1,0))
      }
    }else if(ctype[i] == "E"){
      lvls <- sort(unique(unlist(data[,index])))
      if (!is.null(reflevel[i])){
          index_ref <- which(lvls %in% reflevel[i])
          lvls <- c(lvls[-index_ref], lvls[index_ref])}else{
          lvls <- c(lvls[-1], lvls[1])
        }
      for (j in 1:(length(lvls)-1)){
        newnames <- paste(var, "_level", lvls[j], rest, sep='')
        newdf[,newnames] <- sapply(data[,index], function(x) ifelse(x==lvls[j],1,ifelse(x==lvls[length(lvls)],-1,0)))
      }
    }
  }
  return(newdf)
}
