
#/* The expected value of the canberra location */
#def harm(n):
#  h = 0.0;
#  #print("n to int:\t"); print(n);
#  for i in range(1,int(n+1)):
#    h += 1.0 / float(i);
#  return h;

# harm <- function(n) {
#     return(sum(1/c(1:n)))
# }
# in R version: take a vector as input
harm <- function(n) {
  h  <- 0
  for(i in seq_len(n))
    h <- h+1/i
  return(h)
  # return(sum(1/c(1:n)))
  # return(sum(1/n))
}

#def e_harm( n):
#  return 0.5 * harm(floor(float(n) / 2.0));

e_harm <- function(n) {
  return(0.5 * harm(floor(n / 2.0)))
}

#def o_harm(n):
#  return harm(n) - 0.5 * harm(floor(float(n) / 2.0));

o_harm <- function(n) {
  return(harm(n) - 0.5 * harm(floor(n / 2.0)))
}

a_harm <- function(n){
  return( ifelse(n %% 2, o_harm(n), e_harm(n)))
}

#def a_harm(n):
#  return o_harm(n) if  n % 2 else e_harm(n);

#def canberra(x, y):
#    """Returns the Canberra distance between two P-vectors x and y:
#    sum_i(abs(x_i - y_i) / (abs(x_i) + abs(y_i))).
#    """
#    x_arr = np.ascontiguousarray(x, dtype=np.float)
#    y_arr = np.ascontiguousarray(y, dtype=np.float)

#    if x_arr.shape[0] != y_arr.shape[0]:
#        raise ValueError("x, y: shape mismatch")

#    d = np.sum(np.absolute(np.subtract(x_arr, y_arr))/ np.add(np.absolute(x_arr), np.absolute(y_arr)))

#    return d

canberra <- function(x, y){
    return( sum(abs(x-y)/(abs(x) + abs(y))) )
}

    
#def canberra_location(x, y, k=None):
#    """Returns the Canberra distance between two position lists,
#    `x` and `y`. A position list of length P contains the position 
#    (from 0 to P-1) of P elements. k is the location parameter,
#    if k=None will be set to P.
#    """
#    x_arr = np.ascontiguousarray(x, dtype=np.int)
#    y_arr = np.ascontiguousarray(y, dtype=np.int)

#    if x_arr.shape[0] != y_arr.shape[0]:
#        raise ValueError("x, y: shape mismatch")
#    
#    if k == None:
#        k = x_arr.shape[0]

#    if k <= 0 or k > x_arr.shape[0]:
#        raise ValueError("k must be in [1, %i]" % x_arr.shape[0])   

#        
#    xx = [z+1 if z +1 < k+1 else k+1 for z in x_arr]
#    yy = [z+1 if z +1 < k+1 else k+1 for z in y_arr]

#    xx =  np.array(xx, dtype=float) 
#    yy =  np.array(yy, dtype=float) 

#    d = np.sum(np.divide(np.absolute(np.subtract(xx,yy)), np.add(xx,yy)))

#    return d

canberra_location <- function(x, y, k=NULL){
#    """Returns the Canberra distance between two position lists,
#    `x` and `y`. A position list of length P contains the position 
#    (from 0 to P-1) of P elements. k is the location parameter,
#    if k=None will be set to P.
#    """

    if(length(x) != length(y))
       stop("x, y: shape mismatch")
    
    if(is.null(k))
        k = length(x)

    if(k <= 0 | k > length(x)) {
        stop(paste0("k must be in [1, ", length(x), "\n"))
    }
    xx <- ifelse(x+1 < k+1, x+1, k+1)
    yy <- ifelse(y+1 < k+1, y+1, k+1)    


    d=sum(abs(xx-yy)/(xx+yy))
    return(d)
}

#def canberra_stability(x, k=None):
#    """Returns the Canberra stability indicator between N position
#    lists, where `x` is an (N, P) matrix. A position list of length 
#    P contains the position (from 0 to P-1) of P elements. k is 
#    the location parameter, if k=None will be set to P. The lower 
#    the indicator value, the higher the stability of the lists.

#    The stability is computed by the mean distance of all the 
#    (N(N-1))/2 non trivial values of the distance matrix (computed
#    by canberra_location()) scaled by the expected (average) 
#    value of the Canberra metric.

#    Example:

#    >>> import numpy as np
#    >>> import mlpy
#    >>> x = np.array([[2,4,1,3,0], [3,4,1,2,0], [2,4,3,0,1]])  # 3 position lists
#    >>> mlpy.canberra_stability(x, 3) # stability indicator
#    0.74862979571499755
#    """
#    x_arr = np.ascontiguousarray(x, dtype=np.int)
#    
#    if k == None:
#        k = x_arr.shape[1]

#    if k <= 0 or k > x_arr.shape[1]:
#        raise ValueError("k must be in [1, %i]" % x_arr.shape[1])

#    d = 0.0

#    for i in range(0, x_arr.shape[0]):
#        for j in range(i+1, x_arr.shape[0]):

#            d += canberra_location(x_arr[i], x_arr[j], k) # 2 position lists + location

#    expected = canberra_expected(x_arr.shape[1], k)

#    return (d / ((x_arr.shape[0]*(x_arr.shape[0]-1)) / 2.0) )/expected

canberra_stability <- function(x, k=NULL){
#    """Returns the Canberra stability indicator between N position
#    lists, where `x` is an (N, P) matrix. A position list of length 
#    P contains the position (from 0 to P-1) of P elements. k is 
#    the location parameter, if k=None will be set to P. The lower 
#    the indicator value, the higher the stability of the lists.

#    The stability is computed by the mean distance of all the 
#    (N(N-1))/2 non trivial values of the distance matrix (computed
#    by canberra_location()) scaled by the expected (average) 
#    value of the Canberra metric.

#    Example:

#    >>> import numpy as np
#    >>> import mlpy
#    >>> x = np.array([[2,4,1,3,0], [3,4,1,2,0], [2,4,3,0,1]])  # 3 position lists
#    >>> mlpy.canberra_stability(x, 3) # stability indicator
#    0.74862979571499755
#    """
    
    if(is.null(k))
        k = ncol(x)

    if(k <= 0 | k > ncol(x))
        stop(paste0("k must be in [1, ", ncol(x), "]\n"))

    d <- 0.0

    for(i in 1:(nrow(x)-1)){
        for(j in (i+1):nrow(x)){
            d <- d+canberra_location(x[i,], x[j,], k)
        }
    }

    expected = canberra_expected(ncol(x), k)

    return( (d / ((nrow(x)*(nrow(x)-1)) / 2.0) )/expected )
}
#def canberra_location_expected(p, k=None):
#    """Returns the expected value of the Canberra location distance,
#    where `p` is the number of elements and `k` is the number of 
#    positions to consider.
#    """

#    if k == None:
#        k = p

#    if k <= 0 or k > p:
#        raise ValueError("k must be in [1, %i]" % p)

#    return canberra_expected(p,k)

canberra_location_expected <- function(p, k=NULL){
#    """Returns the expected value of the Canberra location distance,
#    where `p` is the number of elements and `k` is the number of 
#    positions to consider.
#    """
    if(is.null(k)){
        k <- p
    }
    if(k <= 0 | k > p){
        stop(paste0("k must be in [1, %i]",p, "\n"))
    }
    return(canberra_expected(p,k))
}

#def canberra_expected(n, k):
#    mysum = 0.0;
#    for t in range(1,k+1):
#        mysum += t * (a_harm(2 * k - t) - a_harm(t));
#    return 2.0 / n * mysum + (2.0 * (n - k) / n) * (2 * (k + 1) * (harm(2 * k + 1) - harm(k + 1)) - k);
canberra_expected <- function(n, k){
    
    mysum <- 0
    
    for(t in 1:k) {
      mysum <- mysum + t * (a_harm(2 * k - t) - a_harm(t));
    }
    # as implemented, the <>harm functions take only 1 integer !!! cannot vectorize !
    # mysum <- sum( tvect * (a_harm(2 * k - tvect) - a_harm(tvect)) )
    return(2.0 / n * mysum + (2.0 * (n - k) / n) * (2 * (k + 1) * (harm(2 * k + 1) - harm(k + 1)) - k))
}

############################################################################################

# cat("hello\n")

# canberra_expected(5,2);

#x = np.array([[2,4,1,3,0], [3,4,1,2,0], [2,4,3,0,1]])  # 3 position lists
#c_s = canberra_stability(x, 3) # stability indicator
#print(c_s)

#x = matrix(c(2,4,1,3,0,3,4,1,2,0,2,4,3,0,1), byrow = T, ncol=5)

#canberra_stability(x,3)
#canberra_stability(x)

#canberra_expected(1,2)
#import numpy as np
#import mlpy
#x = np.array([[2,4,1,3,0], [3,4,1,2,0], [2,4,3,0,1]])  # 3 position lists
#mlpy.canberra_stability(x, 3) # stability indicator
#0.74862979571499755

