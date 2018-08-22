# create priors from distributions implemented in R
createPriors <- function(pname, shape, p, n, redraw = NULL, fix = NULL){
  priors <- matrix(NA, nrow = n, ncol = length(pname))
  colnames(priors) <- c(pname)
  
  # length checks
  if(length(pname) != length(shape) | 
     length(pname) != length(p) |
     length(shape) != length(p)) stop("Provide as many parameter names, as shape functions, as parameter sets")
  
  for(i in seq_along(pname)){
    if(i %in% fix){
      priors[, pname[i]] <- rep(p[[i]], n)
    }else{
      priors[, pname[i]] <- do.call(shape[[i]], as.list(c(n, p[[i]])))
    }
  }
  
  if(!is.null(redraw)){
    for(i in redraw){
      tofix <- which(priors[,i] < 0)
      fillwith <- priors[tofix, i]
      while(sum(fillwith < 0) > 0){
        newdraw <- do.call(shape[[i]], as.list(c(n, p[[i]])))
        positive <- newdraw[which(newdraw > 0)]
        ifelse(length(positive) > sum(fillwith < 0),
               fillwith[fillwith < 0] <- positive[1:sum(fillwith < 0)],
               fillwith[fillwith < 0][1:length(positive)] <- positive)
      }
      priors[tofix, i] <- fillwith
    }
  }

  return(priors)
}




