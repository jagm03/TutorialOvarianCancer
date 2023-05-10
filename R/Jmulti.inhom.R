### Crossed J function
"Jmulti.inhom" <-
  function(X, I, J, 
           lambda=NULL, lambdaI=NULL, lambdaJ=NULL,
           lambdamin=NULL,
           ...,
           r=NULL, 
           ReferenceMeasureMarkSetI=NULL,
           ratio=FALSE){
    ## compute multitype inhomogeneous G (including determination of r and lmin)
    GIJ <- GmultiInhom(X, I, J,
                       lambda, lambdaI, lambdaJ,
                       lambdamin, ...,
                       r,
                       ReferenceMeasureMarkSetI,
                       ratio)
    ## extract auxiliary values to be used for Finhom
    r <- GIJ$r
    FJ <- FmultiInhom(X, J,
                lambda, lambdaJ, 
                lambdamin,
                ...,
                r=r)
    ## evaluate inhomogeneous J function
    if(!ratio) {
      JIJ <- eval.fv((1 - GIJ) / (1 - FJ))
    } else {
      num <- eval.fv(1 - GIJ)
      den <- eval.fv(1 - FJ)
      JIJ <- eval.fv(num / den)
      JIJ <- rat(JX, num, den)
    }
    ## relabel the fv object
    JIJ <- rebadge.fv(JIJ, quote(J[cross.inhom](r)), c("J","cross.inhom"),
                     names(JIJ), new.labl = attr(GIJ, "labl"))
    ## tack on extra info
    attr(JIJ, "G") <- GIJ
    attr(JIJ, "F") <- FJ
    attr(JIJ, "dangerous") <- attr(GIJ, "dangerous")
    return(JIJ)
  }

"Jdot.inhom" <-
  function(X, i,
           lambda=NULL, lambdaI=NULL, lambdadot=NULL,
           lambdamin=NULL,
           ...,
           r=NULL, 
           ReferenceMeasureMarkSetI = ReferenceMeasureMarkSetI,
           ratio = FALSE){
    X <- as.ppp(X)
    if(!is.marked(X))
      stop(paste("point pattern has no", sQuote("marks")))
    stopifnot(is.multitype(X))
    #
    marx <- marks(X, dfok = FALSE)
    if(missing(i)) i <- levels(marx)[1]
    #  
    I <- (marx == i)
    if(sum(I) == 0)
      stop(paste("No points have mark = ", i))          
    J <- rep.int(TRUE, X$n)
    #  
    result <- Jmulti.inhom(X, I, J, 
                           lambda, lambdaI, lambdadot,
                           lambdamin,
                           ...,
                           r = r, 
                           ReferenceMeasureMarkSetI = ReferenceMeasureMarkSetI,
                           ratio = ratio)
    conserve <- attr(result, "conserve")
    result <- rebadge.as.dotfun(result, "Jinhom", NULL, i)
    attr(result, "conserve") <- conserve
    return(result)
  }

"Jcross.inhom" <- 	
  function(X, i, j, 
           lambda = NULL, lambdaI = NULL, lambdaJ = NULL, 
           lambdamin = NULL,
           ...,
           r = NULL, 
           ReferenceMeasureMarkSetI = NULL,
           ratio = NULL) {
    verifyclass(X, "ppp")
    if(!is.multitype(X, dfok=FALSE))
      stop("Point pattern must be multitype")
    marx <- marks(X)
    if(missing(i))
      i <- levels(marx)[1]
    if(missing(j))
      j <- levels(marx)[2]
    I <- (marx == i)
    J <- (marx == j)
    result <- Jmulti.inhom(X, I, J, 
                           lambda, lambdaI, lambdaJ,
                           lambdamin,
                           ...,
                           r, 
                           ReferenceMeasureMarkSetI,
                           ratio)
    conserve <- attr(result, "conserve")
    result <- rebadge.as.crossfun(result, "Jinhom", NULL, i, j)
    attr(result, "conserve") <- conserve
    return(result)
  }

