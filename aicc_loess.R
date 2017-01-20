#returns a list of aicc and span values as well as the time it took
aicc_opt <- function(in_df, spans = (0.01 + 0.1*(0:9)), resolution = 0.001, cut_factor = 0.1){
  aicc_v <- sapply(spans, aicc, euc_dist = in_df$distance, pos = in_df$pos)
  df <- data.frame(spans, aicc_v)
  current_resolution <- min(diff(spans))
  if (current_resolution > resolution){
    min_vals <- minTwo(localMin(df$aicc_v))
    #gets first column of dataframe after selecting for rows that match the two lowest localMins
    min_spans <- df[df$aicc_v %in% min_vals, 1]
    add_vector <- ((current_resolution * cut_factor) * 1:9)
    new_spans <- c()
    for (x in min_spans){
      new_spans <- append(new_spans, x - add_vector)
      new_spans <- append(new_spans, x + add_vector)
    }
    new_spans <- new_spans[new_spans > 0]
    new_spans <- new_spans[new_spans <= 1]
    #NOTE: will break here if resolution is finer than 0.01 (maybe not?)
    new_spans <- round(new_spans, digits = 3)
    new_spans <- unique(new_spans)
    message("# new spans = ", length(new_spans))
    return(rbind(df,aicc_opt(in_df, spans = new_spans, resolution = resolution, cut_factor = cut_factor)))
  }
  else return(df)
}

#returns minimum two elements (useful for aicc_opt)
minTwo <- function(x){
  len <- length(x)
  if(len<2){
    warning('len < 2; returning x')
    return(x)
  }
  sort(x,partial=c(1, 2))[c(1, 2)]
}

#returns all local minima (problem if repeated local maxima on end)
localMin <- function(x){
  indices <- which(diff(c(FALSE,diff(x)>0,TRUE))>0)
  return(x[indices])
}

aicc <- function (s, euc_dist, pos) {
  # extract values from loess object
  x <- try(loess(euc_dist ~ pos, span=s, degree=1, family="symmetric", surface='direct'), silent=T)
  if(class(x)=="try-error"){return(NA)}
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum( x$residuals^2 ) / (n-1)
  delta1 <- x$one.delta
  delta2 <- x$two.delta
  enp <- x$enp
  #return aicc value
  return(log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2))
  #what I understand:
  #return(n * log(sigma2) + 2*traceL + 2*traceL*(traceL + 1) / (n - traceL - 1))
}