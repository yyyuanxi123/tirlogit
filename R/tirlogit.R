#' Analyze various types of discrete choice data with the Generalized Rank-Ordered Logit Model
#'
#' @description The `tirlogit` function is used to fit the generalized rank-ordered logit model to analyze
#' various types of discrete choice data with and without ties, including the best choice(s),
#' best-worst choice(s), full ranking data with and without ties, etc.
#'
#' The function supports a three-part formula syntax similar to that of the `mlogit` package, which includes:
#' - alternative specific variables with generic coefficients
#' - choice set specific or individual specific variables with specific coefficients
#' - alternative specific variables with alternative specific coefficients.
#'
#' @name tirlogit
#' @import combinat
#' @import roptim
#' @import randtoolbox
#' @import tidyr
#' @param formula a symbolic description of the model to be fitted.
#' @param data an ordinary `data.frame` containing the variables in the model.
#' @param intercept logical; whether the alternative-specific intercepts are included in the model.
#' @param reflevel the reference level (the base alternative) for the utilities when estimating the model
#' (i.e. the one for which the coefficients of choice set specific covariates are set to zero);
#' The default `reflevel` is the alternative whose name appears first when sorted alphabetically from A to Z.
#' @param type a character string specifying the type of the input choice data; "`rank`" for full (complete) rankings,
#' "`best`" for best choice(s), "`BW`" for best-worst choice(s), "`string`" for general rankings in string format.
#' @param shape a character string specifying the shape of the input `data`, either `"long"` or `"wide"`.
#' By default, `shape = "wide"`.
#' @param alts.names only relevant if `shape = "long"`; a character string specifying
#' the name of the column that contains the names of the alternatives.
#' @param indexes a character vector specifying the names of
#' the columns that contain the individual IDs and the choice task IDs.
#' `indexes` can be omitted if `shape = "wide"`.
#' @param ... additional arguments (currently not used).
#' @return An object of class `tirlogit` containing model fit results including the estimated parameters,
#' log-likelihood, the variance-covariance matrix of the estimated parameters, and the convergence status.
#' @export
tirlogit <- function(formula, data, intercept = FALSE, reflevel = NULL,
                     type = "rank", shape = "wide", alts.names = NULL,
                     indexes = NULL,...){
  library(combinat)
  library(randtoolbox)
  library(tidyr)
  
  threshold.dist= NULL
  num.draws = NULL
  ######### functions:
  ## compute the probabilities of complete rankings
  rank_prob=function(t,util,delta){
    M=length(t)
    
    # functions
    seq.logit <- function(v) {
      prod <- apply(v, 1, function(x){
        cum_sums <- rev(cumsum(rev(x)))
        return(prod(x/cum_sums))})
      return(prod)
    }
    
    partial_prob=function(delta, delta.coef, util.shift){
      delta.seq=delta*delta.coef
      util.seq=util.shift+delta.seq[col(util.shift)]
      return(sum(seq.logit(exp(util.seq))))
    }
    
    list <- c()
    for (i in 1:max(t)){
      if(length(which(t %in% i))>1){
        list[[i]] <- permn(which(t %in% i))
      }
      if(length(which(t %in% i))==1){
        list[[i]] <- which(t %in% i)
      }
    }
    #print(list)
    rankgrid=expand.grid(list)
    #colnames(rankgrid)=paste('rank',c(1:max(t))) #c('rank1','rank2','rank3')
    
    rank.convert=matrix(0,nrow=nrow(rankgrid),ncol=M)
    util.shift=matrix(0,nrow=nrow(rankgrid),ncol=M)
    for (i in 1:nrow(rankgrid)){
      rank.convert[i,]=unlist(rankgrid[i,])
      util.shift[i,]=util[,rank.convert[i,]]
    }
    
    diff=diff(sort(as.numeric(t)))
    
    #if there are no ties (a complete ranking)
    if (sum(diff==0)==0){
      delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
      util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
      complete=sum(seq.logit(exp(util.complete)))
      return(complete)
    }
    
    # if there are ties
    if (sum(diff==0)!=0){
      index=which(diff %in% 0)
      delta.coef1=sort(as.numeric(t))-sort(as.numeric(t))[1]
      unprob=partial_prob(delta,delta.coef=delta.coef1,util.shift)
      
      #complete prob without ties
      delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
      util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
      complete=sum(seq.logit(exp(util.complete)))
      
      nposit=length(t)-length(unique(as.numeric(t)))
      
      #t=c(1,1,2,2)
      #t=c(1,1,1,1)
      if(nposit>=2){ # more than one ties
        partprobs=rep(0,(nposit-1))
        for (i in 1:(nposit-1)){# i denotes the number of "?"
          #choose(M-1,i)
          #i=1
          nchoice=choose(nposit,i)
          delmat=matrix(0,ncol=M,nrow=nchoice)
          diff=diff(sort(as.numeric(t)))
          index=which(diff %in% 0)
          combn=combn(index,i)
          
          for (r in 1:nchoice){
            for (c in 2:M){
              delmat[r,c]=delmat[r,c-1]+ifelse((c-1) %in% combn[,r],0,1)
            }
          }
          
          probs=rep(0,nrow(delmat))
          for (k in 1:nrow(delmat)){
            probs[k]=partial_prob(delta,delta.coef=delmat[k,],util.shift)
          }
          
          partprobs[i]=sum(probs)
          
        }
        
        
        tied.prob= unprob-sum(c(rev(partprobs),complete)*rep_len(c(1,-1),nposit))
      }
      if(nposit<2){ # one tie
        tied.prob=unprob-complete}
      return(tied.prob)}
  }
  
  ## compute the probabilities of incomplete rankings
  general_rank_prob=function(t,util,delta,vec){
    
    M=length(t)
    seq.logit <- function(v) {
      prod <- apply(v, 1, function(x){
        cum_sums <- rev(cumsum(rev(x)))
        return(prod(x/cum_sums))})
      return(prod)
    }
    
    partial_prob=function(delta, delta.coef, util.shift){
      delta.seq=delta*delta.coef
      util.seq=util.shift+delta.seq[col(util.shift)]
      return(sum(seq.logit(exp(util.seq))))
    }
    
    if (sum(vec)==0){ # if there are no ? in the ranking
      
      list <- c()
      for (i in 1:max(t)){
        if(length(which(t %in% i))>1){
          list[[i]] <- permn(which(t %in% i))
        }
        if(length(which(t %in% i))==1){
          list[[i]] <- which(t %in% i)
        }
      }
      #print(list)
      rankgrid=expand.grid(list)
      #colnames(rankgrid)=paste('rank',c(1:max(t))) #c('rank1','rank2','rank3')
      
      rank.convert=matrix(0,nrow=nrow(rankgrid),ncol=M)
      util.shift=matrix(0,nrow=nrow(rankgrid),ncol=M)
      for (i in 1:nrow(rankgrid)){
        rank.convert[i,]=unlist(rankgrid[i,])
        util.shift[i,]=util[,rank.convert[i,]]
      }
      
      diff=diff(sort(as.numeric(t)))
      
      #if there are no ties (a complete ranking)
      if (sum(diff==0)==0){
        delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
        util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
        complete=sum(seq.logit(exp(util.complete)))
        return(complete)
      }
      
      # if there are ties
      if (sum(diff==0)!=0){
        index=which(diff %in% 0)
        delta.coef1=sort(as.numeric(t))-sort(as.numeric(t))[1]
        unprob=partial_prob(delta,delta.coef=delta.coef1,util.shift)
        
        #complete prob without ties
        delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
        util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
        complete=sum(seq.logit(exp(util.complete)))
        
        nposit=length(t)-length(unique(as.numeric(t)))
        
        #t=c(1,1,1,2)
        #t=c(2,1,1,1)
        if(nposit>=2){ # more than one ties
          partprobs=rep(0,(nposit-1))
          for (i in 1:(nposit-1)){# i denotes the number of "?"
            #choose(M-1,i)
            nchoice=choose(nposit,i)
            delmat=matrix(0,ncol=M,nrow=nchoice)
            diff=diff(sort(as.numeric(t)))
            index=which(diff %in% 0)
            combn=combn(index,i)
            
            for (r in 1:nchoice){
              for (c in 2:M){
                delmat[r,c]=delmat[r,c-1]+ifelse((c-1) %in% combn[,r],0,1)
              }
            }
            
            probs=rep(0,nrow(delmat))
            for (k in 1:nrow(delmat)){
              probs[k]=partial_prob(delta,delta.coef=delmat[k,],util.shift)
            }
            
            partprobs[i]=sum(probs)
            
          }
          
          
          tied.prob= unprob-sum(c(rev(partprobs),complete)*rep_len(c(1,-1),nposit))
        }
        if(nposit<2){ # one tie
          tied.prob=unprob-complete}
        #return(tied.prob)}
        
      }
    }
    
    
    if (sum(vec)!=0){ # if there exists ? in the ranking
      
      list <- c()
      for (i in 1:max(t)){
        if(length(which(t %in% i))>1){
          list[[i]] <- permn(which(t %in% i))
        }
        if(length(which(t %in% i))==1){
          list[[i]] <- which(t %in% i)
        }
      }
      #print(list)
      rankgrid=expand.grid(list)
      #colnames(rankgrid)=paste('rank',c(1:max(t))) #c('rank1','rank2','rank3')
      
      rank.convert=matrix(0,nrow=nrow(rankgrid),ncol=M)
      util.shift=matrix(0,nrow=nrow(rankgrid),ncol=M)
      for (i in 1:nrow(rankgrid)){
        rank.convert[i,]=unlist(rankgrid[i,])
        util.shift[i,]=util[,rank.convert[i,]]
      }
      
      diff=diff(sort(as.numeric(t)))
      
      #if there are no ties (a complete ranking)
      #if (sum(diff==0)==0){
      #delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
      #util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
      #complete=sum(seq.logit(exp(util.complete)))
      #return(complete)
      #}
      
      # if there are ties (= or ?)
      if (sum(diff==0)!=0){
        #index=which(diff %in% 0)
        delta.coef1=sort(as.numeric(t))-sort(as.numeric(t))[1]
        unprob=partial_prob(delta,delta.coef=delta.coef1,util.shift)
        
        #complete prob without ties
        #delta.complete=delta*(sort(as.numeric(rank.convert[1,]))-sort(as.numeric(rank.convert[1,]))[1])
        #util.complete=util.shift+delta.complete[col(util.shift)] #add a vector to all rows of a matrix
        #complete=sum(seq.logit(exp(util.complete)))
        
        #diff=diff(sort(as.numeric(t)))
        index1=which(diff %in% 0)
        index2=which(vec %in% 1)
        index=index1[!index1 %in% index2]
        nposit=length(index)
        #nposit=length(t)-length(unique(as.numeric(t)))
        
        #t=c(1,1,1,2,2)
        #t=c(1,1,1,2)
        #t=c(2,1,1,1)
        
        
        if(nposit>=2){ # more than one =
          partprobs=rep(0,(nposit-1))
          for (i in 1:(nposit-1)){# i denotes the number of "?"
            #choose(M-1,i)
            nchoice=choose(nposit-1,i)
            delmat=matrix(0,ncol=M,nrow=nchoice)
            #diff=diff(sort(as.numeric(t)))
            #index1=which(diff %in% 0)
            #index2=which(vec %in% 1)
            #index=index[!index1 %in% index2]
            combn=combn(index,i)
            
            #index2=which(vec %in% 1)
            
            for (r in 1:nchoice){
              for (c in 2:M){
                delmat[r,c]=delmat[r,c-1]+ifelse((c-1) %in% combn[,r]|(c-1) %in% index2,0,1)
              }
            }
            
            probs=rep(0,nrow(delmat))
            for (k in 1:nrow(delmat)){
              probs[k]=partial_prob(delta,delta.coef=delmat[k,],util.shift)
            }
            
            partprobs[i]=sum(probs)
            
            #}
            
            
            tied.prob= unprob-sum(rev(partprobs)*rep_len(c(1,-1),nposit-1))
            #return(tied.prob)
          }}
        
        if(nposit==0){
          tied.prob=unprob}
        
        if(nposit==1){
          delmat=rep(0,M)
          for (c in 2:M){
            delmat[c]=delmat[c-1]+ifelse((c-1) %in% index2,0,1)
          }
          tprob=partial_prob(delta,delta.coef=delmat,util.shift)
          tied.prob=unprob-tprob
        }
        #return(tied.prob)}
        #return(tied.prob)
      }}
    
    return(tied.prob)
    
    
  }
  
  
  ###### extract names of variables from the formula
  extractFormula <- function(formula) {
    
    if (!inherits(formula, "formula")) {
      stop("Input must be a formula.")
    }
    
    formula_string <- deparse(formula, width.cutoff = 500)
    
    parts <- strsplit(formula_string, "\\s*~\\s*")[[1]]
    y_part <- parts[1]
    
    rhs_parts <- strsplit(parts[2], "\\s*\\|\\s*")[[1]]
    
    #part <- rhs_parts[3]
    process_part <- function(part) {
      factor_vars <- NULL
      if (part == "0") {
        return(list(variables = NULL, factors = factor_vars))
      } else if (grepl("\\+", part)){
        if (grepl("as.factor", part)){
          matches <- gregexpr("as\\.factor\\(([^)]+)\\)", part, perl = TRUE)
          factor_vars <- regmatches(part, matches)
          factor_vars <- unlist(lapply(factor_vars, function(match) gsub("as\\.factor\\(([^)]+)\\)", "\\1", match)))
          expr_without_factors <- gsub("\\s*\\+?\\s*as\\.factor\\([^)]+\\)\\s*\\+?\\s*", " + ", part)
          expr_without_factors <- gsub("[[:space:]]+|\\+", "", expr_without_factors)
          
          #expr_without_factors <- trimws(gsub("^\\s*\\+\\s*|\\s*\\+\\s*$", "", expr_without_factors))
          #expr_without_factors <- gsub("\\s*\\+\\s*\\+\\s*", " + ", expr_without_factors)
          vars <- unlist(strsplit(expr_without_factors, " \\+ "))
          
          #expr_without_factors <- gsub("as\\.factor\\([^)]+\\) \\+ ", "", part)
          #vars <- unlist(strsplit(expr_without_factors, " \\+ "))
          return(list(variables = vars, factors = factor_vars))
        } else {
          return(list(variables = strsplit(part, "\\s*\\+\\s*")[[1]],  factors = factor_vars))
        }
      } else{
        if (grepl("as.factor", part)){
          matches <- gregexpr("as\\.factor\\(([^)]+)\\)", part, perl = TRUE)
          factor_vars <- regmatches(part, matches)
          factor_vars <- unlist(lapply(factor_vars, function(match) gsub("as\\.factor\\(([^)]+)\\)", "\\1", match)))
          expr_without_factors <- gsub("as\\.factor\\([^)]+\\)( \\+ )?", "", part)
          vars <- unlist(strsplit(expr_without_factors, " \\+ "))
          return(list(variables = vars, factors = factor_vars))
        }else{
          return(list(variables = part,  factors = factor_vars))
        }
        
      }
    }
    
    x_parts <- NULL
    z_parts <- NULL
    w_parts <- NULL
    x_factors <- NULL
    z_factors <- NULL
    w_factors <- NULL
    
    if (length(rhs_parts) >= 1) {
      x_parts <- process_part(rhs_parts[1])$variables
      x_factors <- process_part(rhs_parts[1])$factors
    }
    
    if (length(rhs_parts) >= 2) {
      z_parts <- process_part(rhs_parts[2])$variables
      z_factors <- process_part(rhs_parts[2])$factors
    }
    
    if (length(rhs_parts) == 3) {
      w_parts <- process_part(rhs_parts[3])$variables
      w_factors <- process_part(rhs_parts[3])$factors
    }
    
    return(list(y = y_part, x = x_parts, z = z_parts, w = w_parts,
                x_factors = x_factors, z_factors = z_factors, w_factors = w_factors))
    
  }
  
  
  extractRank <- function(rank_string, names_alternatives) {
    
    rank_string <- gsub(" ", "", rank_string)
    
    letters <- unlist(strsplit(rank_string, "[^A-Za-z]"))
    letters <- letters[letters != ""]
    
    if (!setequal(letters, names_alternatives)) {
      stop("The ranking strings are not specified correctly.")
    }
    
    symbols <- unlist(strsplit(rank_string, "[A-Za-z]"))
    symbols <- symbols[symbols != ""]
    
    equal_signs <- sum(symbols == '=') ## the number of equal signs "="
    vec <- as.numeric(symbols == "?")
    rank_changes <- as.numeric(symbols == ">")
    ranking <- cumsum(c(1, rank_changes))
    names(ranking) <- letters
    ranking <- ranking[names_alternatives]
    
    list(ranking = ranking, vec = vec, equal_signs = equal_signs)
  }
  
  
  #### a function to extract data
  extractData <- function(formula, data, reflevel){
    
    ## names of alternatives, reference level and ranking matrix
    
    if (type != "string") {
      
      if (type == 'rank'){
        rank_columns <- grep(paste("^", formula_parts$y, "\\.", sep=""), names(data), value = TRUE)
        names_alternatives <- sub(paste("^", formula_parts$y, "\\.", sep=""), "", rank_columns)
        names_alternatives <- sort(names_alternatives) # Alphabetical order A-Z
        
        if(is.null(reflevel)) {reflevel <- names_alternatives[1]}  # default reflevel
        index_ref <- which(names_alternatives %in% reflevel)  # the index of the reflevel
        names_alternatives <- c(names_alternatives[-index_ref], reflevel)
        
        ## ranking matrix
        rank_names <- paste(formula_parts$y, names_alternatives, sep="." )
        rank_matrix <- as.matrix(data[, rank_names])
        vec_matrix <- NULL
        
        ## ties detection
        has_ties <- apply(rank_matrix, 1, function(row) any(duplicated(row)))
      }
      
      if (type == 'best'){
        rank_columns <- grep(paste("^", formula_parts$y, "\\.", sep=""), names(data), value = TRUE)
        names_alternatives <- sub(paste("^", formula_parts$y, "\\.", sep=""), "", rank_columns)
        names_alternatives <- sort(names_alternatives)
        
        if(is.null(reflevel)) {reflevel <- names_alternatives[1]}
        #if(is.null(reflevel)) {reflevel <- names_alternatives[length(names_alternatives)]}  # default reflevel
        index_ref <- which(names_alternatives %in% reflevel)  # the index of the reflevel
        names_alternatives <- c(names_alternatives[-index_ref], reflevel)
        
        ## ranking matrix
        rank_names <- paste(formula_parts$y, names_alternatives, sep="." )
        rank_matrix <- as.matrix(data[, rank_names])
        num_best <- rowSums(rank_matrix)
        
        ## ties detection
        has_ties <- num_best > 1
        
        # convert rank matrix
        rank_matrix[rank_matrix==0] <- 2
        vec_matrix <- matrix(1, nrow = nrow(rank_matrix), ncol = ncol(rank_matrix)-1)
        num_quesMark <- ncol(rank_matrix) - num_best - 1
        num_quesMark[num_quesMark<0] <- 0
        for (i in 1: nrow(rank_matrix)){
          indicator <- ncol(rank_matrix) - num_quesMark[i]- 1
          vec_matrix[i, 1:indicator] <- 0
        }
      }
    }
    
    
    if (type == 'BW'){
      rank_columns <- grep(paste("^", formula_parts$y, "\\.", sep=""), names(data), value = TRUE)
      names_alternatives <- sub(paste("^", formula_parts$y, "\\.", sep=""), "", rank_columns)
      names_alternatives <- sort(names_alternatives)
      
      if(is.null(reflevel)) {reflevel <- names_alternatives[1]}
      index_ref <- which(names_alternatives %in% reflevel)  # the index of the reflevel
      names_alternatives <- c(names_alternatives[-index_ref], reflevel)
      num_alts <- length(names_alternatives)
      
      ## ranking matrix
      rank_names <- paste(formula_parts$y, names_alternatives, sep="." )
      rank_matrix <- as.matrix(data[, rank_names])
      num_best <- rowSums(rank_matrix == 1)
      num_worst <- rowSums(rank_matrix == -1)
      num_mid <- num_alts - num_best - num_worst
      
      ## ties detection
      has_ties <- num_best > 1 | num_worst > 1
      
      
      rank_matrix[rank_matrix == 0] <- 2
      rank_matrix[rank_matrix == -1] <- 3
      for (i in which(num_mid == 0)) {
        rank_matrix[i, rank_matrix[i, ] == 3] <- 2
      }
      
      vec_matrix <- matrix(0, nrow = nrow(rank_matrix), ncol = num_alts-1)
      num_quesMark <- num_alts - num_best - num_worst - 1
      num_quesMark[num_quesMark<0] <- 0
      for (i in 1:nrow(rank_matrix)){
        if(num_quesMark[i] > 0){
          vec_matrix[i, (num_best[i]+1):(num_best[i]+num_quesMark[i])] <- 1
        }
        #indicator <- num_alts - num_quesMark[i] - 1
        #vec_matrix[i, 1:indicator] <- 0
      }
    }
    
    
    if (type == "string") {
      rank_extract <- data[1,formula_parts$y]
      rank_extract <- gsub(" ", "", rank_extract)
      names_alternatives <- unlist(strsplit(rank_extract, "[^A-Za-z]"))
      names_alternatives <- names_alternatives[names_alternatives != ""]
      names_alternatives <- sort(names_alternatives)
      
      if(is.null(reflevel)) {reflevel <- names_alternatives[1]}
      #if(is.null(reflevel)) {reflevel <- names_alternatives[length(names_alternatives)]}  # default reflevel
      index_ref <- which(names_alternatives %in% reflevel)  # the index of the reflevel
      names_alternatives <- c(names_alternatives[-index_ref], reflevel)
      
      result_list <- lapply(data[,formula_parts$y], function(rank_string) extractRank(rank_string, names_alternatives = names_alternatives))
      rank_matrix <- do.call(rbind, lapply(result_list, function(item) item$ranking))
      vec_matrix <- do.call(rbind, lapply(result_list, function(item) item$vec))
      
      ## ties detection
      has_ties <- as.vector(do.call(rbind, lapply(result_list, function(item) item$equal_signs))) > 0
    }
    
    ## alternative vars. with generic coefficients
    vars_alt_generic <- formula_parts$x
    factors_alt_generic <- formula_parts$x_factors
    
    x_matrices_list_org <- list()
    x_factors_list_org <- list()
    
    
    if (is.null(vars_alt_generic) | length(vars_alt_generic) == 0){x_matrices_list_org = NULL
    }else{
      for (i in 1: length(vars_alt_generic)){
        xi_columns <- paste(vars_alt_generic[i], names_alternatives, sep="." )
        x_matrices_list_org[[i]] <- as.matrix(data[, xi_columns])
      }
    }
    
    
    if (is.null(factors_alt_generic) | length(factors_alt_generic) == 0){
      x_factors_list_org <- NULL
      names_x_factors <- NULL
    }else{
      names_x_factors <- c()
      
      for (k in 1: length(factors_alt_generic)){
        xi_columns <- paste(factors_alt_generic[k], names_alternatives, sep="." )
        newdf <- data.frame(matrix(nrow = nrow(data), ncol = 0))
        lvls <- sort(unique(unlist(data[,xi_columns])))
        names <- paste0("as.factor(",factors_alt_generic[k],")", lvls[1:(length(lvls)-1)])
        names_x_factors <- c(names, names_x_factors)
        
        for (j in 1:(length(lvls)-1)){
          newnames <- paste(paste0("as.factor(",factors_alt_generic[k],")", lvls[j]),  names_alternatives, sep = '.')
          newdf[,newnames] <- sapply(data[,xi_columns], function(x) ifelse(x==lvls[j],1,0))
        }
        x_factors_list_org[[k]] <- newdf
        #names_x_factors <- unlist(names_x_factors)
      }
    }
    #names_x_factors[2]
    
    ## individual vars. with specific coefficients
    vars_indv_specific <- formula_parts$z
    factors_indv_specific <- formula_parts$z_factors
    data_indv_specific_org <- as.matrix(data[,vars_indv_specific])
    z_factors_list_org <- list()
    
    if (is.null(factors_indv_specific)){
      z_factors_list_org <- NULL
      names_z_factors <- NULL
    }else{
      names_z_factors <- c()
      for (k in 1: length(factors_indv_specific)){
        newdf <- data.frame(matrix(nrow = nrow(data), ncol = 0))
        lvls <- sort(unique(unlist(data[,factors_indv_specific[k]])))
        names <- paste0("as.factor(",factors_indv_specific[k],")", lvls[1:(length(lvls)-1)])
        names_z_factors <- c(names, names_z_factors)
        
        for (j in 1:(length(lvls)-1)){
          newnames <- paste0("as.factor(",factors_indv_specific[k],")", lvls[j])
          newdf[,newnames] <- sapply(data[,factors_indv_specific[k]], function(x) ifelse(x==lvls[j],1,0))
        }
        z_factors_list_org[[k]] <- newdf
        #names_z_factors <- unlist(names_z_factors)
      }
    }
    
    
    
    ## alternative vars. with specific coefficients
    vars_alt_specific <- formula_parts$w
    factors_alt_specific <- formula_parts$w_factors
    
    w_matrices_list_org <- list()
    w_factors_list_org <- list()
    
    if(is.null(vars_alt_specific) | length(vars_alt_specific) == 0) {
      w_matrices_list_org <- NULL
    }else{
      for (i in 1: length(vars_alt_specific)){
        wi_columns <- paste(vars_alt_specific[i], names_alternatives, sep="." )
        w_matrices_list_org[[i]] <- as.matrix(data[, wi_columns])
        #w_matrices_list[[i]] <- (w_matrices_list_org[[i]] - mean(unlist(w_matrices_list_org[[i]])))/sd(unlist(w_matrices_list_org[[i]]))
      }}
    
    if (is.null(factors_alt_specific)){
      w_factors_list_org = NULL
      names_w_factors <- NULL
    }else{
      names_w_factors <- c()
      for (k in 1: length(factors_alt_specific)){
        wi_columns <- paste(factors_alt_specific[k], names_alternatives, sep="." )
        newdf <- data.frame(matrix(nrow = nrow(data), ncol = 0))
        lvls <- sort(unique(unlist(data[,wi_columns])))
        names <- paste0("as.factor(",factors_alt_specific[k],")", lvls[1:(length(lvls)-1)])
        names_w_factors <- c(names_w_factors, names)
        
        for (j in 1:(length(lvls)-1)){
          #names_w_factors[[k]][j] <-  paste0("as.factor(",factors_alt_specific[k],")", lvls[j])
          newnames <- paste(paste0("as.factor(",factors_alt_specific[k],")", lvls[j]),  names_alternatives, sep = '.')
          newdf[,newnames] <- sapply(data[,wi_columns], function(x) ifelse(x==lvls[j],1,0))
        }
        w_factors_list_org[[k]] <- newdf
        #names_w_factors <- unlist(names_w_factors)
      }
    }
    
    
    return(list(names_alternatives = names_alternatives,
                rank_matrix = rank_matrix,
                vec_matrix = vec_matrix,
                x_matrices_list_org = x_matrices_list_org,
                x_factors_list_org = x_factors_list_org,
                data_indv_specific_org = data_indv_specific_org,
                z_factors_list_org = z_factors_list_org,
                w_matrices_list_org = w_matrices_list_org,
                w_factors_list_org = w_factors_list_org,
                names_x_factors = names_x_factors,
                names_z_factors = names_z_factors,
                names_w_factors = names_w_factors,
                reflevel = reflevel,
                has_ties = has_ties))
  }
  
  
  
  long2wide <- function(formula_parts, data, indexes, alt){
    
    #formula_parts <- extractFormula(formula)
    
    df_wide <- pivot_wider(data, id_cols = all_of(c(indexes, formula_parts$z)),
                           names_from = alt,
                           values_from = all_of(c(formula_parts$x, formula_parts$w, formula_parts$y)),
                           names_sep = ".")
    df_wide <- as.data.frame(df_wide)
    
    return(df_wide)
  }
  
  
  ####### main ######
  formula_parts <- extractFormula(formula)
  
  #is.null(all_of(formula_parts == NULL))
  if(shape == 'long'){
    data <- long2wide(formula_parts, data, indexes, alt.names)
  }
  
  num_task <- nrow(data)
  extract_data <- extractData(formula, data, reflevel)
  reflevel <- extract_data$reflevel
  names_alternatives <- extract_data$names_alternatives
  rank_matrix <- extract_data$rank_matrix
  vec_matrix <- extract_data$vec_matrix
  x_matrices_list_org <- extract_data$x_matrices_list_org
  x_factors_list_org <- extract_data$x_factors_list_org
  data_indv_specific_org <- extract_data$data_indv_specific_org
  z_factors_list_org <- extract_data$z_factors_list_org
  w_matrices_list_org <- extract_data$w_matrices_list_org
  w_factors_list_org <- extract_data$w_factors_list_org
  names_x_factors <- extract_data$names_x_factors
  ### factors: not yet delete the reference group!!!!
  names_z_factors <- extract_data$names_z_factors
  names_w_factors <- extract_data$names_w_factors
  
  has_ties <- extract_data$has_ties
  ties <- (sum(has_ties) > 0)
  
  num_alt <- length(names_alternatives)
  
  num_alt_generic <- length(formula_parts$x)
  num_indv_specific <- length(formula_parts$z)
  num_alt_specific <- length(formula_parts$w)
  
  num_alt_generic_factors <- length(formula_parts$x_factors)
  #length(names_x_factors)
  num_indv_specific_factors <- length(names_z_factors)
  #length(formula_parts$z_factors)
  #length(names_z_factors)
  num_alt_specific_factors <- length(formula_parts$w_factors)
  #length(names_w_factors)
  #length(formula_parts$w_factors)
  #length(names_w_factors)
  
  
  
  
  
  x_matrices_list <- NULL
  w_matrices_list <- NULL
  data_indv_specific <- data_indv_specific_org
  #x_factors_list <- NULL
  # w_factors_list <- NULL
  # z_factors_list <- NULL
  
  start_idx_1 <- 1
  end_idx_1 <- length(x_matrices_list_org)
  
  if (length(x_factors_list_org) == 0){
    start_idx_factor_1 <- end_idx_1 + 1
    end_idx_factor_1 <- end_idx_1
  }else{
    start_idx_factor_1 <- end_idx_1 + 1
    end_idx_factor_1 <- end_idx_1 + sum(sapply(x_factors_list_org, ncol)) / num_alt
    #end_idx_factor_1 <- end_idx_1 + length(x_factors_list_org) * sum(sapply(x_factors_list_org, ncol)) / num_alt
  }
  
  start_idx_2 <- end_idx_factor_1 + 1
  end_idx_2 <- end_idx_factor_1 + (num_alt-1) * ncol(data_indv_specific_org)
  
  if (length(z_factors_list_org) == 0){
    start_idx_factor_2 <- end_idx_2 + 1
    end_idx_factor_2 <- end_idx_2
  }else{
    start_idx_factor_2 <- end_idx_2 + 1
    end_idx_factor_2 <- end_idx_2 + (num_alt-1) * sum(sapply(z_factors_list_org, ncol))
  }
  
  start_idx_3 <- end_idx_factor_2 + 1
  end_idx_3 <- end_idx_factor_2 + num_alt * length(w_matrices_list_org)
  if (length(z_factors_list_org) == 0){
    start_idx_factor_3 <- end_idx_3 + 1
    end_idx_factor_3 <- end_idx_3
  }else{
    start_idx_factor_3 <- end_idx_3 + 1
    end_idx_factor_3 <- end_idx_3 + sum(sapply(w_factors_list_org, ncol))
  }
  
  
  
  if (intercept == TRUE){
    
    num_para <- (num_alt-1) + end_idx_factor_3
    start_idx_4 <- end_idx_factor_3 + 1
    end_idx_4 <- num_para
    names_para <- rep(NA, num_para)
    mean_std <- rep(0, num_para)
    sd_std <- rep(1, num_para)
    
    #change_intercept_indv <- matrix(0, nrow = end_idx_factor_2-start_idx_2+1, ncol = num_alt-1)
    #change_intercept_alt <- matrix(0, nrow = end_idx_factor_3-start_idx_3+1, ncol = num_alt-1)
    
    change_intercept_indv <- matrix(0, nrow = num_indv_specific, ncol = num_alt-1)
    change_intercept_alt <- matrix(0, nrow = num_alt_specific, ncol = num_alt-1)
    
    
    ### standardization
    if (start_idx_1 <= end_idx_1){
      x_matrices_list <- lapply(x_matrices_list_org, function(matrix) {
        (matrix - mean(unlist(matrix))) / sd(unlist(matrix))
      })
      
      mean_std[start_idx_1:end_idx_1] <- sapply(x_matrices_list_org, function(x) mean(unlist(x)))
      sd_std[start_idx_1:end_idx_1] <-  sapply(x_matrices_list_org, function(x) sd(unlist(x)))
      
      ### names of the estimated parameters
      names_para[start_idx_1:end_idx_1] <- formula_parts$x
      
    }
    
    
    if (start_idx_factor_1 <= end_idx_factor_1){
      
      #x_factors_list <- lapply(x_factors_list_org, function(matrix) {
      # (matrix - mean(unlist(matrix))) / sd(unlist(matrix))
      #  })
      
      # mean_std[start_idx_factor_1:end_idx_factor_1] <- sapply(x_factors_list_org, function(x) mean(unlist(x)))
      # sd_std[start_idx_factor_1:end_idx_factor_1] <- sapply(x_factors_list_org, function(x) sd(unlist(x)))
      
      ### names of the estimated parameters
      names_para[start_idx_factor_1:end_idx_factor_1] <- names_x_factors
      
    }
    
    #change_intercept_indv_list <- list()
    if (start_idx_2 <= end_idx_2){
      data_indv_specific <- apply(data_indv_specific_org, MARGIN = 2, function(x) {
        (x - mean(x)) / sd(x)
      })
      #num_indv_1 <- ncol(data_indv_specific_org)
      mean_std[start_idx_2:end_idx_2] <- rep(apply(data_indv_specific_org, MARGIN = 2, mean), each = num_alt - 1)
      sd_std[start_idx_2:end_idx_2] <- rep(apply(data_indv_specific_org, MARGIN = 2, sd), each = num_alt - 1)
      
      change_intercept_indv <- matrix(mean_std[start_idx_2:end_idx_2]/sd_std[start_idx_2:end_idx_2],
                                      nrow = num_indv_specific, byrow = TRUE)
      names_para[start_idx_2:end_idx_2] <- unlist(sapply(1:length(formula_parts$z), function(i) {
        paste(formula_parts$z[i], names_alternatives[-num_alt], sep=':')
      }))
      #change_intercept_indv_list[[1]] <- matrix(mean_std[start_idx_2:end_idx_2]/sd_std[start_idx_2:end_idx_2],
      #  nrow = num_indv_1, ncol = num_alt-1, byrow = TRUE)
      
      # names_para[start_idx_2:end_idx_2] <- unlist(sapply(1:num_indv_1, function(i) {
      #  paste(formula_parts$z[i], names_alternatives[-num_alt], sep=':')
      # }))
    }
    
    #colnames(z_factors_list_org[[1]])
    if (start_idx_factor_2 <= end_idx_factor_2){
      
      #z_factors_list <- lapply(z_factors_list_org, function(matrix) {
      #(matrix - mean(unlist(matrix))) / sd(unlist(matrix))
      #})
      
      #z_factors_list <- apply(z_factors_list_org, MARGIN = 2, function(x) {
      #(x - mean(x)) / sd(x)
      #})
      
      # mean_std[start_idx_factor_2:end_idx_factor_2] <- rep(sapply(z_factors_list_org, function(x) mean(unlist(x))),  each = num_alt - 1)
      #rep(apply(data_indv_specific_org, MARGIN = 2, mean), each = num_alt - 1)
      #  sd_std[start_idx_factor_2:end_idx_factor_2] <- rep(sapply(z_factors_list_org, function(x) sd(unlist(x))),  each = num_alt - 1)
      #rep(apply(data_indv_specific_org, MARGIN = 2, sd), each = num_alt - 1)
      # change_intercept_indv_list[[2]] <- matrix(mean_std[start_idx_factor_2:end_idx_factor_2]/sd_std[start_idx_factor_2:end_idx_factor_2],
      #       nrow = length(names_z_factors) , ncol = num_alt-1, byrow = TRUE)
      
      names_para[start_idx_factor_2:end_idx_factor_2] <- unlist(sapply(1:length(names_z_factors), function(i) {
        paste(names_z_factors[i], names_alternatives[-num_alt], sep=':')
      }))
    }
    
    
    #change_intercept_indv <- do.call(rbind, change_intercept_indv_list)
    
    if (start_idx_3 <= end_idx_3) {
      
      w_matrices_list <- lapply(w_matrices_list_org, function(matrix) {
        (matrix - mean(unlist(matrix))) / sd(unlist(matrix))
      })
      
      mean_std[start_idx_3:end_idx_3] <- rep(sapply(w_matrices_list_org, function(x) mean(unlist(x))), each = num_alt)
      sd_std[start_idx_3:end_idx_3] <-  rep(sapply(w_matrices_list_org, function(x) sd(unlist(x))), each = num_alt)
      change_intercept_alt <- matrix(rep(sapply(w_matrices_list_org, function(x) mean(unlist(x)/sd(unlist(x)))), each = num_alt -1),
                                     nrow = length(formula_parts$w), byrow = TRUE)
      names_para[start_idx_3:end_idx_3] <- unlist(sapply(1:length(formula_parts$w), function(i) {
        paste(formula_parts$w[i], names_alternatives, sep=':')
      }))
    }
    
    
    if (start_idx_factor_3 <= end_idx_factor_3){
      names_para[start_idx_factor_3:end_idx_factor_3] <- unlist(sapply(1:length(names_w_factors), function(i) {
        paste(names_w_factors[i], names_alternatives, sep=':')
      }))
    }
    
    
    names_para[start_idx_4:end_idx_4] <- paste('(Intercept)', names_alternatives[-num_alt], sep=':')
    
    
  }else{
    
    num_para <- end_idx_factor_3
    names_para <- rep(NA, num_para)
    sd_std <- rep(1, num_para)
    
    
    ### standardization
    if (start_idx_1 <= end_idx_1){
      
      x_matrices_list <- lapply(x_matrices_list_org, function(matrix) {
        matrix / sd(unlist(matrix))
      })
      sd_std[start_idx_1:end_idx_1] <-  sapply(x_matrices_list_org, function(x) sd(unlist(x)))
      names_para[start_idx_1:end_idx_1] <- formula_parts$x
      
    }
    
    
    if (start_idx_factor_1 <= end_idx_factor_1){
      
      names_para[start_idx_factor_1:end_idx_factor_1] <- names_x_factors
      
    }
    
    
    if (start_idx_2 <= end_idx_2){
      
      data_indv_specific <- apply(data_indv_specific_org, MARGIN = 2, function(x) {
        x / sd(x)
      })
      sd_std[start_idx_2:end_idx_2] <- rep(apply(data_indv_specific_org, MARGIN = 2, sd), each = num_alt - 1)
      names_para[start_idx_2:end_idx_2] <- unlist(sapply(1:num_indv_specific, function(i) {
        paste(formula_parts$z[i], names_alternatives[-num_alt], sep=':')
      }))
      
    }
    
    if (start_idx_factor_2 <= end_idx_factor_2){
      
      names_para[start_idx_factor_2:end_idx_factor_2] <- unlist(sapply(1:length(names_z_factors), function(i) {
        paste(names_z_factors[i], names_alternatives[-num_alt], sep=':')
      }))
    }
    
    
    if (start_idx_3 <= end_idx_3) {
      w_matrices_list <- lapply(w_matrices_list_org, function(matrix) {
        matrix / sd(unlist(matrix))
      })
      
      sd_std[start_idx_3:end_idx_3] <-  rep(sapply(w_matrices_list_org, function(x) sd(unlist(x))), each = num_alt)
      
      names_para[start_idx_3:end_idx_3] <- unlist(sapply(1:num_alt_specific, function(i) {
        paste(formula_parts$w[i], names_alternatives, sep=':')
      }))
    }
    
    if (start_idx_factor_3 <= end_idx_factor_3){
      names_para[start_idx_factor_3:end_idx_factor_3] <- unlist(sapply(1:length(names_w_factors), function(i) {
        paste(names_w_factors[i], names_alternatives, sep=':')
      }))
    }
    
    
  }
  
  
  #util_array <- array(0, dim = c(num_task, num_alt, (num_alt_generic + num_indv_specific + num_alt_specific +
  #                                      ifelse(num_alt_generic_factors>0,1,0) + num_indv_specific_factors + ifelse(num_alt_specific_factors>0,1,0) )))
  
  list_data <- list(NULL)
  
  if((num_alt_generic + num_alt_generic_factors) > 0){
    list_data[[1]] <- do.call(cbind,c(x_matrices_list, x_factors_list_org))
  } else{
    list_data[[1]] <- NULL
  }
  #list_data[[1]] <- do.call(cbind,c(NULL, NULL))
  
  
  
  
  zpara_length <- num_indv_specific + num_indv_specific_factors
  zref_index <- seq(num_alt, by = num_alt, length.out = zpara_length)
  zorg_index <- 1:(zpara_length * num_alt)
  zpara_index <- setdiff(zorg_index, zref_index)
  
  if(length(z_factors_list_org) > 0) {
    data_indv_factors <- do.call(cbind, z_factors_list_org)
  } else{
    data_indv_factors <- NULL
  }
  
  #data_indv_factors <- do.call(cbind, z_factors_list_org)
  zdata <- cbind(data_indv_specific, data_indv_factors)
  if((num_indv_specific + num_indv_specific_factors) > 0){
    list_data[[2]] <- zdata[,rep(1:zpara_length, each = num_alt)]
  } else{
    list_data[[2]] <- NULL
  }
  
  
  if((num_alt_specific + num_alt_specific_factors) >0){
    list_data[[3]] <- do.call(cbind,c(w_matrices_list, w_factors_list_org))
  } else{
    list_data[[3]] <- NULL
  }
  
  
  #rep(1:ypara_length, each = num_alt)
  #combined_data <- c(xdata,zdata)
  #head(xdata)
  #head(zdata)
  combined_data <- do.call(cbind,list_data)
  
  
  
  
  #cbind(xdata,zdata,wdata)
  #head(combined_data)
  #combined_data <- do.call(cbind,append(xdata, zdata, wdata))
  #ncol(combined_data)
  
  #xpara <- rep(para[start_idx_1:end_idx_factor_1], each = num_alt)
  #zpara <- rep(0, zpara_length * num_alt)
  #zpara[zpara_index] <- para[start_idx_2:end_idx_factor_2]
  #wpara <- para[start_idx_3:end_idx_factor_3]
  
  #combined_para <- c(xpara, zpara, wpara)
  
  
  #y_length <- num_indv_specific + num_indv_specific_factors
  #ypara_length <- (num_indv_specific + num_indv_specific_factors) * num_alt
  #rep(0, ypara_length)
  #seq(1, ypara_length, by = 3)
  #para[start_idx_2:end_idx_factor_2]
  
  # t(t(combined_data) * combined_para)[1:5,]
  #xpara * xdata[1,]
  
  
  if (ties == FALSE){
    ### compute loglikeliood -84.0388
    #para <- m1$par  -84.0388 < -84.03878
    #para <- ml.MC3$coefficients[names_para]*sd_std
    #m2$par/sd_std
    
    f <- function(para){
      
      if(start_idx_1 <= end_idx_factor_1){
        xpara <- rep(para[start_idx_1:end_idx_factor_1], each = num_alt)
      } else {
        xpara <- NULL
      }
      if(start_idx_2 <= end_idx_factor_2){
        zpara <- rep(0, zpara_length * num_alt)
        zpara[zpara_index] <- para[start_idx_2:end_idx_factor_2]
      } else {
        zpara <- NULL
      }
      if(start_idx_3 <= end_idx_factor_3){
        wpara <- para[start_idx_3:end_idx_factor_3]
      } else{
        wpara <- NULL
      }
      combined_para <- c(xpara, zpara, wpara)
      
      combined_util <- t(t(combined_data) * combined_para)
      
      lengthout <- ncol(combined_util)/num_alt
      util_matrix <- matrix(0, nrow = num_task, ncol = num_alt)
      for (i in 1:num_alt){
        seq_index <- seq(i, by = num_alt, length.out = lengthout)
        if (length(seq_index) == 1){util_matrix[,i] = combined_util[, seq_index]} else{
          util_matrix[,i] = rowSums(combined_util[, seq_index])}
        #util_matrix[,i] = rowSums(combined_util[, seq_index])
      }
      
      #head(util_matrix)
      #newmat <- lapply(util_matrix, function(i) {util_matrix[,i] = rowSums(combined_util[, seq(i, by = num_alt, length.out = lengthout)])})
      #newmat
      
      ## intercept
      if (intercept == TRUE){
        coef_intercept <- para[start_idx_4:end_idx_4]
        coef_intercept <- c(coef_intercept, 0)
        util_intercept <- rep(1,nrow(data))%o% coef_intercept
        util_matrix <- util_matrix + util_intercept
        #apply(util_array, c(1, 2), sum) + util_intercept
      }
      
      ### end for !is.null(combinded_para)
      #if(intercept == FALSE){
      #util_sum <-
      #apply(util_array, c(1, 2), sum)
      #}
      
      probIndv <- rep(NA, num_task)
      if (is.null(vec_matrix)){
        for (i in 1:num_task){
          probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0)
        }
      }else{
        for (i in 1:num_task){
          probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0, vec = vec_matrix[i,])
        }
      }
      
      probIndv[abs(probIndv) < 1e-10] <- 1e-10
      f <- -sum(log(probIndv),na.rm = TRUE)
    }
    
    
    f_null <- function(para){
      #coef_intercept <- para[start_idx_4:end_idx_4]
      coef_intercept <- c(para, 0)
      util_matrix <- rep(1,nrow(data))%o% coef_intercept
      #util_matrix <- util_intercept
      
      probIndv <- rep(NA, num_task)
      if (is.null(vec_matrix)){
        for (i in 1:num_task){
          probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0)
        }
      }else{
        for (i in 1:num_task){
          probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0, vec = vec_matrix[i,])
        }
      }
      
      probIndv[abs(probIndv) < 1e-10] <- 1e-10
      f <- -sum(log(probIndv),na.rm = TRUE)
    }
    
    # para <- para.o
    para.o <- rep(0,num_para) #rnorm(num_para, 0, 1)
    m1 <- optim(para.o,f,hessian=TRUE,method="L-BFGS-B", control = list(maxit = 20000))
    
    
    
    #### compute fitted.values ########
    estimate_m1 <- m1$par
    combined_util <- t(t(combined_data) * estimate_m1[start_idx_1:end_idx_3])
    lengthout <- ncol(combined_util)/num_alt
    util_matrix <- matrix(0, nrow = num_task, ncol = num_alt)
    for (i in 1:num_alt){
      seq_index <- seq(i, by = num_alt, length.out = lengthout)
      if (length(seq_index) == 1){util_matrix[,i] = combined_util[, seq_index]} else{
        util_matrix[,i] = rowSums(combined_util[, seq_index])}
    }
    if (intercept == TRUE){
      coef_intercept <- estimate_m1[start_idx_4:end_idx_4]
      coef_intercept <- c(coef_intercept, 0)
      util_intercept <- rep(1,nrow(data))%o% coef_intercept
      util_matrix <- util_matrix + util_intercept
    }
    
    probIndv <- rep(NA, num_task)
    if (is.null(vec_matrix)){
      for (i in 1:num_task){
        probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0)
      }
    }else{
      for (i in 1:num_task){
        probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=0, vec = vec_matrix[i,])
      }
    }
    
    fitted.values <- probIndv
    
    if (intercept == TRUE){
      mat_intercept <- matrix(0, nrow = num_para, ncol = num_alt-1)
      diag(mat_intercept[start_idx_4: end_idx_4, ]) <- rep(1, (num_alt-1))
      
      
      if (num_indv_specific > 0){
        para_z <- matrix(m1$par[start_idx_2:end_idx_2], byrow= TRUE, nrow = num_indv_specific)
        #minus_intercept[start_idx_4:end_idx_4] <- colSums(para_z * change_intercept_indv)
        
        for (i in 1:num_indv_specific){
          diag(mat_intercept[(start_idx_2 + (i-1)*(num_alt-1)):(start_idx_2 + i*(num_alt-1) -1), ]) <-
            - change_intercept_indv[i,]
        }
      }
      
      if (num_alt_specific > 0){
        para_w <- matrix(m1$par[start_idx_3:end_idx_3], byrow= TRUE, nrow = num_alt_specific)
        para_w <- para_w - para_w[,num_alt]
        para_w <- para_w[,-num_alt]
        #minus_intercept[start_idx_4:end_idx_4] <- minus_intercept[start_idx_4:end_idx_4] +
        #colSums(change_intercept_alt * para_w)
        
        for (i in 1:num_indv_specific){
          diag(mat_intercept[(start_idx_3 + (i-1)*(num_alt-1)):(start_idx_3 + i*(num_alt-1) -1), ]) <-
            - change_intercept_alt[i,]
          mat_intercept[start_idx_3 + i*(num_alt-1) , ] <- change_intercept_alt[i,]
        }
        
      }
      
      intercept_org <- colSums(m1$par*mat_intercept)
      
      vcov <- solve(m1$hessian)
      
      vcov_converted <- vcov
      sd_stds <- c(sd_std,1)
      for (i in 1:end_idx_3){
        for (j in 1:end_idx_3){
          vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
        }
      }
      vcov_converted_1<-vcov_converted[start_idx_1:end_idx_3, start_idx_1:end_idx_3]
      
      #vcov_converted_1<-vcov_converted[start_idx_1:end_idx_3,start_idx_1:end_idx_3]
      if(end_idx_3 - start_idx_1 >=1){
        colnames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
        rownames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]}else{
          names(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
        }
      
      
      #sqrt(diag(vcov_converted))
      
      vcov_converted_2<-t(mat_intercept)%*%vcov%*%mat_intercept
      colnames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
      rownames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
      
      estimate <- c(m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3], intercept_org)
      stdError <- c(sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3], diag(sqrt(t(mat_intercept)%*%vcov%*%mat_intercept)))
      z_value <- estimate/stdError # H0:coef=0
      # two-tailed p-value
      p_value=rep(0.5,length(z_value))
      for (i in 1:length(z_value)){
        if (z_value[i]>=0){
          p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
        }
        else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
        }
      }
      
      #p_value <- format.pval(object$result$p_value, digits = 5, eps = .Machine$double.eps)
      #p_value <- round(p_value,6)
      
      result <- cbind.data.frame(Estimate=estimate,StdError=stdError,z_value, p_value)
      #result <- cbind.data.frame(Estimate=(m1$par - minus_intercept)/sd_std,StdError=sqrt(diag(vcov))/sd_std,z_value, p_value)
      rownames(result) <- names_para
      names(estimate) <- names_para
      
      
      
      para.null <- rep(0,end_idx_4-start_idx_4 +1) #rnorm(num_para, 0, 1)
      m_null <- optim(para.null,f_null,hessian=TRUE,method="L-BFGS-B", control = list(maxit = 20000))
      
      LL_null = -m_null$value
      LL_model = -m1$value
      chisq_stat <- 2 * (LL_model - LL_null)
      p_chisqTest <- pchisq(chisq_stat, df = end_idx_3, lower.tail = FALSE)
      McFadden_R2 <- 1-LL_model/LL_null
      AIC_value <- -2 * LL_model + 2 * length(estimate)
      BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
      output <- list(call = match.call(),
                     result = result,
                     coefficients = estimate,
                     LogLikelihood = LL_model,
                     vcov = list(vcov_converted_1, vcov_converted_2),
                     LL_null = LL_null,
                     McFadden_R2 = McFadden_R2,
                     chisq_stat = chisq_stat,
                     p_chisqTest = p_chisqTest,
                     fitted.values = fitted.values,
                     AIC_value = AIC_value,
                     BIC_value = BIC_value,
                     ties = ties,
                     convergence = (m1$convergence==0))
    }else{
      vcov <- solve(m1$hessian)
      
      vcov_converted <- vcov
      sd_stds <- c(sd_std,1)
      for (i in 1:nrow(vcov_converted)){
        for (j in 1:nrow(vcov_converted)){
          vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
        }
      }
      
      colnames(vcov_converted) <- names_para
      rownames(vcov_converted) <- names_para
      
      estimate <- m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3]
      stdError <- sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3]
      z_value <- estimate/stdError # H0:coef=0
      # two-tailed p-value
      p_value=rep(0.5,length(z_value))
      for (i in 1:length(z_value)){
        if (z_value[i]>=0){
          p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
        }
        else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
        }
      }
      
      #p_value <- round(p_value,6)
      result <- cbind.data.frame(Estimate=estimate,StdError=stdError,z_value, p_value)
      rownames(result) <- names_para
      names(estimate) <- names_para
      
      LL_model <- -m1$value
      AIC_value <- -2 * LL_model + 2 * length(estimate)
      BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
      
      output <- list(call = match.call(),
                     result = result,
                     coefficients = estimate,
                     LogLikelihood = -m1$value,
                     vcov = vcov_converted,
                     fitted.values = fitted.values,
                     AIC_value = AIC_value,
                     BIC_value = BIC_value,
                     ties = ties,
                     convergence = (m1$convergence==0))
    }
  }
  
  
  if(!is.null(threshold.dist)){
    #id <- data
    task_mat <-  data[,indexes]
    colnames(task_mat) <- c("id","taskid")
    #task_mat$id <- data[,indexes[1]]
    #task_mat$taskid <- data[,indexes[2]]
    id <- unique(task_mat[,1])
    
    #threshold.dist =NULL
    
    if(threshold.dist == 'lognorm'){
      #library(randtoolbox)
      R = num.draws
      halton_seq <- halton(n = R, dim = 1)
      norm_data <- qnorm(halton_seq)
      #LikeliMat <- matrix(0, nrow = I, ncol = R)
      #delta_mean <- exp()
      #para <- c(model_health3$coefficients[1:5], -2, 0.2)
      #delta <- exp(para[num_para+1] + para[num_para+2] * norm_data)
      #delta <- rep(log(0.12259014), R)
      #para = c(estimate, 0.2)
      
      #delta <- rep(0.5,10)
      f <- function(para){
        
        delta <- exp(para[num_para+1] + para[num_para+2] * norm_data)
        
        if(start_idx_1 <= end_idx_factor_1){
          xpara <- rep(para[start_idx_1:end_idx_factor_1], each = num_alt)
        } else {
          xpara <- NULL
        }
        if(start_idx_2 <= end_idx_factor_2){
          zpara <- rep(0, zpara_length * num_alt)
          zpara[zpara_index] <- para[start_idx_2:end_idx_factor_2]
        } else {
          zpara <- NULL
        }
        if(start_idx_3 <= end_idx_factor_3){
          wpara <- para[start_idx_3:end_idx_factor_3]
        } else{
          wpara <- NULL
        }
        combined_para <- c(xpara, zpara, wpara)
        
        combined_util <- t(t(combined_data) * combined_para)
        
        lengthout <- ncol(combined_util)/num_alt
        util_matrix <- matrix(0, nrow = num_task, ncol = num_alt)
        for (i in 1:num_alt){
          seq_index <- seq(i, by = num_alt, length.out = lengthout)
          if (length(seq_index) == 1){util_matrix[,i] = combined_util[, seq_index]} else{
            util_matrix[,i] = rowSums(combined_util[, seq_index])}
        }
        
        ## intercept
        if (intercept == TRUE){
          coef_intercept <- para[start_idx_4:end_idx_4]
          coef_intercept <- c(coef_intercept, 0)
          util_intercept <- rep(1,nrow(data))%o% coef_intercept
          util_matrix <- util_matrix + util_intercept
          #apply(util_array, c(1, 2), sum) + util_intercept
        }
        
        #prob_product <- rep(0, length(id))
        LikeliMat <- matrix(0, nrow = length(id), ncol = R)
        
        for (indv in 1:length(id)){
          
          for (r in 1:R){
            task_rows <- which(task_mat$id == id[indv])
            probIndv <- rep(NA, num_task)
            
            if (is.null(vec_matrix)){
              for (i in task_rows){
                probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=delta[r])
              }
            }else{
              for (i in task_rows){
                probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=delta[r], vec = vec_matrix[i,])
              }
            }
            #probIndv[abs(probIndv) < 1e-10] <- 1e-10
            LikeliMat[indv, r] <- prod(probIndv[task_rows])
          }
        }
        #LikeliMat[LikeliMat < 1e-10]  <- 1e-10
        f <- -sum(log(rowSums(LikeliMat)/R))
        return(f)
      }
      # para <- para.o
      para.o <- c(rep(0,(num_para)),0, 1)
      lower.bounds <- rep(-Inf, (num_para + 2))
      upper.bounds <- rep(Inf,(num_para + 2))
      
      m1 <- optim(para.o,f, hessian=TRUE, method = "L-BFGS-B", lower = lower.bounds, upper = upper.bounds, control = list(maxit = 20000))
      
      if (intercept == TRUE){
        mat_intercept <- matrix(0, nrow = num_para, ncol = num_alt-1)
        diag(mat_intercept[start_idx_4: end_idx_4, ]) <- rep(1, (num_alt-1))
        
        
        if (num_indv_specific > 0){
          para_z <- matrix(m1$par[start_idx_2:end_idx_2], byrow= TRUE, nrow = num_indv_specific)
          #minus_intercept[start_idx_4:end_idx_4] <- colSums(para_z * change_intercept_indv)
          
          for (i in 1:num_indv_specific){
            diag(mat_intercept[(start_idx_2 + (i-1)*(num_alt-1)):(start_idx_2 + i*(num_alt-1) -1), ]) <-
              - change_intercept_indv[i,]
          }
        }
        
        if (num_alt_specific > 0){
          para_w <- matrix(m1$par[start_idx_3:end_idx_3], byrow= TRUE, nrow = num_alt_specific)
          para_w <- para_w - para_w[,num_alt]
          para_w <- para_w[,-num_alt]
          #minus_intercept[start_idx_4:end_idx_4] <- minus_intercept[start_idx_4:end_idx_4] +
          #colSums(change_intercept_alt * para_w)
          
          for (i in 1:num_indv_specific){
            diag(mat_intercept[(start_idx_3 + (i-1)*(num_alt-1)):(start_idx_3 + i*(num_alt-1) -1), ]) <-
              - change_intercept_alt[i,]
            mat_intercept[start_idx_3 + i*(num_alt-1) , ] <- change_intercept_alt[i,]
          }
          
        }
        
        intercept_org <- colSums(m1$par[1:num_para]*mat_intercept)
        
        vcov <- solve(m1$hessian)
        
        vcov_converted <- vcov
        sd_stds <- c(sd_std,1)
        #sd_stds[start_idx_4:end_idx_4] <- diag((t(mat_intercept)%*%vcov%*%mat_intercept))
        for (i in 1:end_idx_3){
          for (j in 1:end_idx_3){
            vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
          }
        }
        vcov_converted_1<-vcov_converted[start_idx_1:end_idx_3, start_idx_1:end_idx_3]
        #vcov_converted_1<-vcov_converted[c(start_idx_1:end_idx_3, num_para+1),c(start_idx_1:end_idx_3, num_para+1)]
        if(end_idx_3-start_idx_1 >= 1){
          colnames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
          rownames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
        }else{
          names(vcov_converted_1)<- names_para[start_idx_1:end_idx_3]
        }
        #sqrt(diag(vcov_converted))
        
        vcov_converted_2<-t(mat_intercept)%*%vcov[1:num_para,1:num_para]%*%mat_intercept
        colnames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
        rownames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
        
        estimate <- c(m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                      intercept_org , m1$par[(num_para+1):(num_para+2)])
        
        stdError <- c(sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                      diag(sqrt(t(mat_intercept)%*%vcov[1:num_para,1:num_para]%*%mat_intercept)), sqrt(diag(vcov))[(num_para+1):(num_para+2)])
        
        
        z_value <- estimate/stdError # H0:coef=0
        # two-tailed p-value
        p_value=rep(0.5,length(z_value))
        for (i in 1:length(z_value)){
          if (z_value[i]>=0){
            p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
          }
          else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
          }
        }
        
        #p_value <- round(p_value,6)
        result <- cbind.data.frame(Estimate=estimate,
                                   StdError=stdError,
                                   z_value, p_value)
        
        rownames(result) <- c(names_para, "mu","sigma")
        names(estimate) <- c(names_para, "mu","sigma")
        
        LL_model <- -m1$value
        AIC_value <- -2 * LL_model + 2 * length(estimate)
        BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
        
        output <- list(call = match.call(),
                       result = result,
                       coefficients = estimate,
                       LogLikelihood = -m1$value,
                       vcov = list(vcov_converted_1, vcov_converted_2),
                       AIC_value = AIC_value,
                       BIC_value = BIC_value,
                       ties = ties,
                       convergence = (m1$convergence==0)
        )
      }else{
        vcov <- solve(m1$hessian)
        
        vcov_converted <- vcov
        sd_stds <- c(sd_std,1)
        for (i in 1:nrow(vcov_converted)){
          for (j in 1:nrow(vcov_converted)){
            vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
          }
        }
        
        colnames(vcov_converted) <- c(names_para, "mu","sigma")
        rownames(vcov_converted) <- c(names_para, "mu","sigma")
        
        
        estimate <- c(m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                      m1$par[(num_para+1):(num_para+2)])
        
        stdError <- c(sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                      sqrt(diag(vcov))[(num_para+1):(num_para+2)])
        
        
        z_value <- estimate/stdError # H0:coef=0
        # two-tailed p-value
        p_value=rep(0.5,length(z_value))
        for (i in 1:length(z_value)){
          if (z_value[i]>=0){
            p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
          }
          else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
          }
        }
        
        #p_value <- round(p_value,6)
        result <- cbind.data.frame(Estimate=estimate,
                                   StdError=stdError,
                                   z_value, p_value)
        
        rownames(result) <- c(names_para, "mu","sigma")
        names(estimate) <- c(names_para, "mu","sigma")
        
        LL_model <- -m1$value
        AIC_value <- -2 * LL_model + 2 * length(estimate)
        BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
        
        output <- list(call = match.call(),
                       result = result,
                       coefficients = estimate,
                       LogLikelihood = -m1$value,
                       vcov = vcov_converted,
                       AIC_value = AIC_value,
                       BIC_value = BIC_value,
                       ties = ties,
                       convergence = (m1$convergence==0)
        )
      }
      
    }
  }
  
  
  #### fixed threshold
  
  if (ties == TRUE & is.null(threshold.dist)){
    
    
    ### compute loglikeliood
    f <- function(para){
      
      
      if(start_idx_1 <= end_idx_factor_1){
        xpara <- rep(para[start_idx_1:end_idx_factor_1], each = num_alt)
      } else {
        xpara <- NULL
      }
      if(start_idx_2 <= end_idx_factor_2){
        zpara <- rep(0, zpara_length * num_alt)
        zpara[zpara_index] <- para[start_idx_2:end_idx_factor_2]
      } else {
        zpara <- NULL
      }
      if(start_idx_3 <= end_idx_factor_3){
        wpara <- para[start_idx_3:end_idx_factor_3]
      } else{
        wpara <- NULL
      }
      combined_para <- c(xpara, zpara, wpara)
      
      combined_util <- t(t(combined_data) * combined_para)
      
      lengthout <- ncol(combined_util)/num_alt
      util_matrix <- matrix(0, nrow = num_task, ncol = num_alt)
      for (i in 1:num_alt){
        seq_index <- seq(i, by = num_alt, length.out = lengthout)
        if (length(seq_index) == 1){util_matrix[,i] = combined_util[, seq_index]} else{
          util_matrix[,i] = rowSums(combined_util[, seq_index])}
      }
      
      ## intercept
      if (intercept == TRUE){
        coef_intercept <- para[start_idx_4:end_idx_4]
        coef_intercept <- c(coef_intercept, 0)
        util_intercept <- rep(1,nrow(data))%o% coef_intercept
        util_matrix <- util_matrix + util_intercept
        #apply(util_array, c(1, 2), sum) + util_intercept
      }
      
      probIndv <- rep(NA, num_task)
      if (is.null(vec_matrix)){
        for (i in 1:num_task){
          probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=para[num_para+1])
        }
      }else{
        for (i in 1:num_task){
          probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=para[num_para+1], vec = vec_matrix[i,])
        }
      }
      
      probIndv[abs(probIndv) < 1e-10] <- 1e-10
      f <- -sum(log(probIndv),na.rm = TRUE)
    }
    
    # para <- para.o
    para.o <- c(rep(0,(num_para)),1)
    lower.bounds <- c(rep(-Inf, num_para), 0)
    upper.bounds <- rep(Inf,(num_para + 1))
    
    #m1 <- optim(para.o,f,hessian=TRUE,method = "L-BFGS-B", lower = lower.bounds, upper = upper.bounds)
    m1 <- optim(para.o,f, hessian=TRUE, method = "L-BFGS-B", lower = lower.bounds, upper = upper.bounds, control = list(maxit = 20000))
    
    #### compute fitted.values ########
    estimate_m1 <- m1$par
    combined_util <- t(t(combined_data) * estimate_m1[start_idx_1:end_idx_3])
    lengthout <- ncol(combined_util)/num_alt
    util_matrix <- matrix(0, nrow = num_task, ncol = num_alt)
    for (i in 1:num_alt){
      seq_index <- seq(i, by = num_alt, length.out = lengthout)
      if (length(seq_index) == 1){util_matrix[,i] = combined_util[, seq_index]} else{
        util_matrix[,i] = rowSums(combined_util[, seq_index])}
    }
    if (intercept == TRUE){
      coef_intercept <- estimate_m1[start_idx_4:end_idx_4]
      coef_intercept <- c(coef_intercept, 0)
      util_intercept <- rep(1,nrow(data))%o% coef_intercept
      util_matrix <- util_matrix + util_intercept
    }
    
    probIndv <- rep(NA, num_task)
    if (is.null(vec_matrix)){
      for (i in 1:num_task){
        probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=estimate_m1[num_para+1])
      }
    }else{
      for (i in 1:num_task){
        probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=estimate_m1[num_para+1], vec = vec_matrix[i,])
      }
    }
    
    fitted.values <- probIndv
    
    if (intercept == TRUE){
      mat_intercept <- matrix(0, nrow = num_para, ncol = num_alt-1)
      diag(mat_intercept[start_idx_4: end_idx_4, ]) <- rep(1, (num_alt-1))
      
      
      if (num_indv_specific > 0){
        para_z <- matrix(m1$par[start_idx_2:end_idx_2], byrow= TRUE, nrow = num_indv_specific)
        #minus_intercept[start_idx_4:end_idx_4] <- colSums(para_z * change_intercept_indv)
        
        for (i in 1:num_indv_specific){
          diag(mat_intercept[(start_idx_2 + (i-1)*(num_alt-1)):(start_idx_2 + i*(num_alt-1) -1), ]) <-
            - change_intercept_indv[i,]
        }
      }
      
      if (num_alt_specific > 0){
        para_w <- matrix(m1$par[start_idx_3:end_idx_3], byrow= TRUE, nrow = num_alt_specific)
        para_w <- para_w - para_w[,num_alt]
        para_w <- para_w[,-num_alt]
        #minus_intercept[start_idx_4:end_idx_4] <- minus_intercept[start_idx_4:end_idx_4] +
        #colSums(change_intercept_alt * para_w)
        
        for (i in 1:num_alt_specific){
          diag(mat_intercept[(start_idx_3 + (i-1)*(num_alt-1)):(start_idx_3 + i*(num_alt-1) -1), ]) <-
            - change_intercept_alt[i,]
          mat_intercept[start_idx_3 + i*(num_alt-1) , ] <- change_intercept_alt[i,]
        }
        
      }
      
      intercept_org <- colSums(m1$par[1:num_para]*mat_intercept)
      
      vcov <- solve(m1$hessian)
      
      vcov_converted <- vcov
      sd_stds <- c(sd_std,1)
      #sd_stds[start_idx_4:end_idx_4] <- diag((t(mat_intercept)%*%vcov%*%mat_intercept))
      for (i in 1:end_idx_3){
        for (j in 1:end_idx_3){
          vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
        }
      }
      vcov_converted_1<-vcov_converted[start_idx_1:end_idx_3, start_idx_1:end_idx_3]
      
      if(end_idx_3-start_idx_1 >= 1){
        colnames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
        rownames(vcov_converted_1) <- names_para[start_idx_1:end_idx_3]
      }else{
        names(vcov_converted_1)<- names_para[start_idx_1:end_idx_3]
      }
      
      vcov_converted_2<-t(mat_intercept)%*%vcov[1:num_para,1:num_para]%*%mat_intercept
      colnames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
      rownames(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
      #names(vcov_converted_2) <- names_para[start_idx_4:length(names_para)]
      #sqrt(diag(vcov_converted_1))
      #sqrt(diag(vcov_converted_2))
      #vcov_converted <- vcov
      #sd_stds <- c(sd_std,1)
      #for (i in 1:nrow(vcov_converted)){
      # for (j in 1:nrow(vcov_converted)){
      #  vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
      # }
      #}
      
      #colnames(vcov_converted) <- c(names_para, "threshold")
      #rownames(vcov_converted) <- c(names_para, "threshold")
      
      estimate <- c(m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                    intercept_org , m1$par[num_para+1])
      
      stdError <- c(sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                    diag(sqrt(t(mat_intercept)%*%vcov[1:num_para,1:num_para]%*%mat_intercept)), sqrt(diag(vcov))[num_para+1])
      
      
      z_value <- estimate/stdError # H0:coef=0
      # two-tailed p-value
      p_value=rep(0.5,length(z_value))
      for (i in 1:length(z_value)){
        if (z_value[i]>=0){
          p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
        }
        else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
        }
      }
      
      #p_value <- round(p_value,6)
      result <- cbind.data.frame(Estimate=estimate,
                                 StdError=stdError,
                                 z_value, p_value)
      
      rownames(result) <- c(names_para, "threshold")
      names(estimate) <- c(names_para, "threshold")
      
      f_null <- function(para){
        
        coef_intercept <- para[1:(end_idx_4-start_idx_4+1)]
        coef_intercept <- c(coef_intercept, 0)
        util_matrix <- rep(1,nrow(data))%o% coef_intercept
        
        #apply(util_array, c(1, 2), sum) + util_intercept
        
        probIndv <- rep(NA, num_task)
        if (is.null(vec_matrix)){
          for (i in 1:num_task){
            probIndv[i] <- rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=para[length(para)])
          }
        }else{
          for (i in 1:num_task){
            probIndv[i] <- general_rank_prob(rank_matrix[i,],util=matrix(util_matrix[i,],nrow=1),delta=para[length(para)], vec = vec_matrix[i,])
          }
        }
        
        probIndv[abs(probIndv) < 1e-10] <- 1e-10
        f <- -sum(log(probIndv),na.rm = TRUE)
      }
      
      para.null <- c(rep(0,end_idx_4-start_idx_4+1), 1) #rnorm(num_para, 0, 1)
      m_null <- optim(para.null,f_null,hessian=FALSE,method="L-BFGS-B", control = list(maxit = 20000))
      
      LL_null = -m_null$value
      LL_model = -m1$value
      chisq_stat <- 2 * (LL_model - LL_null)
      p_chisqTest <- pchisq(chisq_stat, df = end_idx_3, lower.tail = FALSE)
      McFadden_R2 <- 1-LL_model/LL_null
      AIC_value <- -2 * LL_model + 2 * length(estimate)
      BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
      
      
      output <- list(call = match.call(),
                     result = result,
                     coefficients = estimate,
                     LogLikelihood = -m1$value,
                     vcov = list(vcov_converted_1, vcov_converted_2),
                     LL_null = LL_null,
                     McFadden_R2 = McFadden_R2,
                     chisq_stat = chisq_stat,
                     p_chisqTest = p_chisqTest,
                     fitted.values = fitted.values,
                     AIC_value = AIC_value,
                     BIC_value = BIC_value,
                     ties = ties,
                     convergence = (m1$convergence==0))
    }else{
      vcov <- solve(m1$hessian)
      
      vcov_converted <- vcov
      sd_stds <- c(sd_std,1)
      for (i in 1:nrow(vcov_converted)){
        for (j in 1:nrow(vcov_converted)){
          vcov_converted[i,j] <- vcov_converted[i,j]/(sd_stds[i] * sd_stds[j])
        }
      }
      
      colnames(vcov_converted) <- c(names_para, "threshold")
      rownames(vcov_converted) <- c(names_para, "threshold")
      
      
      estimate <- c(m1$par[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                    m1$par[num_para+1])
      
      stdError <- c(sqrt(diag(vcov))[start_idx_1:end_idx_3]/sd_std[start_idx_1:end_idx_3],
                    sqrt(diag(vcov))[num_para+1])
      
      
      z_value <- estimate/stdError # H0:coef=0
      # two-tailed p-value
      p_value=rep(0.5,length(z_value))
      for (i in 1:length(z_value)){
        if (z_value[i]>=0){
          p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = FALSE)
        }
        else{p_value[i]=2*pnorm(z_value[i], mean = 0, sd = 1, lower.tail = TRUE)
        }
      }
      
      #p_value <- round(p_value,6)
      result <- cbind.data.frame(Estimate=estimate,
                                 StdError=stdError,
                                 z_value, p_value)
      
      rownames(result) <- c(names_para, "threshold")
      names(estimate) <- c(names_para, "threshold")
      
      LL_model <- -m1$value
      AIC_value <- -2 * LL_model + 2 * length(estimate)
      BIC_value <- -2 * LL_model + log(nrow(data)) * length(estimate)
      
      output <- list(call = match.call(),
                     result = result,
                     coefficients = estimate,
                     LogLikelihood = -m1$value,
                     vcov = vcov_converted,
                     fitted.values = fitted.values,
                     AIC_value = AIC_value,
                     BIC_value = BIC_value,
                     ties = ties,
                     convergence = (m1$convergence==0))
    }
    
  }
  
  class(output) <- "tirlogit"
  return(output)
}

#' Summary method for tirlogit model objects.
#'
#' @description Provides a summary of the estimation results of the fitted `tirlogit` model.
#' @param object a fitted model of class `tirlogit`.
#' @param ... additional arguments (currently not used).
#' @return Estimation results of the fitted tirlogit model, including:
#' - The model call, showing how the `tirlogit` model was specified.
#' - A table that includes the estimates of the coefficients, and the corresponding standard errors, z-values, p-values,
#'   and significance codes.
#' - The log-likelihood of the fitted model.
#' - McFadden R^2 (only relevant if `Intercept = TRUE`).
#' - Likelihood ratio test of the fitted model (vs. null model, only relevant if `Intercept = TRUE`).
#' @export
summary.tirlogit <- function(object, ...) {
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  if (object$ties == TRUE) {
    cat("Ties detected in the data.\n")
  } else {
    cat("No ties detected in the data; the threshold was fixed at 0.\n")
  }
  cat("\n")
  cat("Coefficients:\n")
  signif_codes <- cut(object$result$p_value,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1),
                      labels = c("***", "**", "*", ".", " "))
  object$result$Signif <- signif_codes
  object$result$p_value <- format.pval(object$result$p_value, digits = 4, eps = .Machine$double.eps)
  object$result$Estimate <- formatC(object$result$Estimate, format = "f", digits = 6)
  object$result$StdError <- formatC(object$result$StdError, format = "f", digits = 6)
  object$result$z_value <- formatC(object$result$z_value, format = "f", digits = 4)
  #object$result$p_value <- formatC(object$result$p_value, format = "e", digits = 4)
  
  colnames(object$result) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)", "")
  print(object$result, row.names = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("\nLog-Likelihood: ", object$LogLikelihood, "\n")
  if (!is.null(object$McFadden_R2)) {
    cat("McFadden R^2:", format(object$McFadden_R2, digits = 4), "\n")
  }
  
  if (!is.null(object$chisq_stat)) {
    cat("Likelihood ratio test: chisq =", format(object$chisq_stat, digits = 5),
        paste0("(p-value = ", formatC(object$p_chisqTest, format = "e", digits = 5), ")\n"))
  }
  invisible(object)
}

#' Information Criteria for the tirlogit model
#'
#' @description Calculates the Akaike Information Criterion (AIC)
#' and Bayesian Information Criterion (BIC) for an object of class `tirlogit`,
#' based on the log-likelihood and the number of parameters.
#' @param object a fitted model of class `tirlogit`.
#' @param ... additional arguments (currently not used).
#' @return For `AIC.tirlogit`, the AIC value. For `BIC.tirlogit`, the BIC value.
#' @examples
#' data("mode_wide", package = "tirlogit")
#' model_modeChoices<- tirlogit(multiBest ~ cost | dist + income , data = mode_wide,
#'                             intercept = FALSE, reflevel = "car",
#'                             ties = TRUE, type = "best", shape = "wide",
#'                             indexes = c("case", "taskid"), alts.names = "alt")
#' summary(model_modeChoices)
#' AIC(model_modeChoices)
#' BIC(model_modeChoices)
#' @name InformationCriteria
#' @export
AIC.tirlogit <- function(object, ...) {
  aic_value <- object$AIC_value
  return(aic_value) }
#' @rdname InformationCriteria
#' @export
BIC.tirlogit <- function(object, ...) {
  bic_value <- object$BIC_value
  return(bic_value)
}

#' Compute Willingness-to-Pay (WTP) Estimates for a tirlogit Model

#'
#' @description Computes the willingness-to-pay (WTP) estimates for specified variables
#' relative to a chosen cost variable in a fitted `tirlogit` model.
#' @param model a fitted model of class `tirlogit`.
#' @param var.names a character vector specifying one or more variable names for
#' which the WTP will be computed.
#' @param wrt a single character string specifying the name of the cost variable
#' (the variable with respect to which the WTP is calculated).
#' @return A table summarizing the WTP estimates, including the estimate, standard error,
#' t-statistic, p-value, and significance codes for each specified variable.
#' @examples
#' mode_long$time <- mode_long$ivt + mode_long$ovt
#' model0 <- tirlogit(formula = best ~ cost + freq + time, data = mode_long, intercept = FALSE,
#'                    reflevel = NULL, type = "best", shape = "long",
#'                    indexes = c("case", "taskid"),
#'                    alts.names = "alt")
#' summary(model0)
#' wtp(model0, var.names = "freq" , wrt = "cost") # Compute WTP for a single variable
#' wtp(model0, var.names = c("freq", "time") , wrt = "cost") # Compute WTP for multiple variables
#' @export
wtp <- function(model, var.names, wrt){
  
  if (!all(c(var.names, wrt) %in% names(model$coefficients))) {
    stop("One or more of the variables (", paste(var.names, collapse = ", "), ", ", wrt, ") are not present in the model coefficients.")
  }
  
  results_list <- list()
  
  for (i in 1:length(var.names)){
    var1 = var.names[i]
    var2 = wrt
    if(is.list(model$vcov)){vcov <- model$vcov[[1]]}else{vcov <- model$vcov}
    beta1 <- model$coefficients[var1]
    beta2 <- model$coefficients[var2]
    var_beta1 <- vcov[var1,var1]
    var_beta2 <- vcov[var2,var2]
    cov_beta1_beta2 <- vcov[var1,var2]
    
    mrs_est <- beta1/beta2
    mrs_var <- (var_beta1 / beta2^2) +
      ((beta1^2 * var_beta2) / beta2^4) -
      (2 * beta1 * cov_beta1_beta2 / beta2^3)
    mrs_se <- sqrt(mrs_var)
    t_value <- mrs_est / mrs_se
    p_value <- 2 * pt(-abs(t_value), df = Inf)
    
    results_list[[i]] <- data.frame(
      Estimate = formatC(mrs_est, format = "f", digits = 6),
      StdError = formatC(mrs_se, format = "f", digits = 6),
      t_value = formatC(t_value, format = "f", digits = 4),
      p_value = format.pval(p_value, digits = 4, eps = .Machine$double.eps),
      row.names = paste(var1, "/", var2, sep = "")
    )
    
  }
  
  result <- do.call(rbind, results_list)
  signif_codes <- cut(p_value,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1),
                      labels = c("***", "**", "*", ".", " "))
  result$Signif <- signif_codes
  colnames(result) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)", "")
  cat("Willigness-to-pay Estimate:\n")
  print(result, row.names = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

