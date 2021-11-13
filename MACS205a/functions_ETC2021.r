dividif=function(x,y){
## Computes the divided differences (coefficients on the Newton basis) for
## Lagrange interpolation.
##
## @title dividif: Newton's Divided differences
## @param x a vector containing the interpolation nodes
## @param y a vector of same size as x: values of the interpolated function at
## the nodes
## @return a vector of same size as x: the divided differences
## \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'.
n = length(x) -1 ## n: degree of Lagrange polynomial.
Tmat = matrix(ncol = n+1, nrow = n+1)
Tmat[,1] = y ## initialisation of the vector of divided differences:
if(n ==0) {return(diag(Tmat))}
for (j in 2:(n+1) ) {
Tmat[j : (n+1), j ] = (Tmat[j : (n+1), j-1 ] - Tmat[(j-1): n, j-1 ])/(x[j : (n+1)] - x[1:(n+2-j)])
}
return(diag(Tmat))
}


hornerNewton = function(a,x,z){
## Horner's method: Evaluates a polynom P at points z, given
## nodes x and the coefficients a of P in Newton's basis
##
## @param a : vector: the coefficients of the polynomial in
## Newton's basis
## @param x : the interpolation nodes.
## @param z : vector of points where the polynom needs to be
## evaluated.
## @return : a vector of same size as z: the value of the
## polynomial at points z.
##
n <- length(x) - 1 ## degree of the Lagrange poynomial
if( (n < 0) || (length(a) != (n+1)) )
{
stop('at least one interpolating point is needed,
a and x should have same length')
}
f <- a[n+1]*rep(1,length(z))
if(n >= 1){
for( k in 1:n){
f <- a[n-k+1] + f*(z-x[n-k+1])
}
}
return(f)
}


interpolDividif=function(x,y,z){
## Efficient Lagrange interpolation using Horner's method with
## Newton basis for evaluation
## @param x : vector containing the interpolation nodes
## @param y : vector of same size as x: values of the interpolated
## function at the nodes
## @param z : vector of points where the interpolating polynomial
## needs to be evaluated.
## @return : vector of same size as z: the value of the
## interpolating polynomial at points z.
a <- dividif(x,y)
return (hornerNewton(a,x,z))
}

tchebychevNodes = function(n) {
  return (cos(((0:(n-1) + 0.5)*pi/n)))
}



interpolLagrange =function(n, a, b, neval, nodes = 'equi', FUN, Plot, Return = TRUE){
## Generic Lagrange interpolation, with equidistant or Chebyshev nodes.
## @param n : the degree of the interpolating polynomial on each
## subinterval
## @param a : left end-point of the interval
## @param b : right end-point of the interval
## @param neval :number of evaluation points (a regular grid will be
## used on [a,b]
## @param nodes :string, either "equi" (default) for equidistant
## Lagrange interpolation (on each subinterval) or "cheby" for
## using Chebyshev nodes.
## @param FUN: the function to be interpolated
## @param Plot : logical. Setting 'Plot' to TRUE produces a plot
## showing the graph of
## the true functions and its interpolation.
## @return : vector of size 'neval': the values of the Lagrange
## polynomial on an equi-distant grid.
if (nodes == "equi"){
x = seq(-1, 1, length.out = n+1)
}
else if (nodes == "cheby"){
x = tchebychevNodes(n+1)
}
else{stop("the nodes must be either 'equi' or 'cheby'") }
##
z = seq(a, b, length.out = neval)
u = (z - (a+b)/2)*2/(b-a)
y = sapply((a+b)/2 + (b-a)/2*x, FUN = FUN)
f = interpolDividif(x, y, u)
##
##
if( Plot ){
    if (nodes == "equi"){ methodName = " equidistant "}
    else { methodName = " Chebyshev "}
    
    plot(z, sapply(z,FUN), type="l", ylim=range(c(y,f)) )
    title(main = paste("Lagrange interpolation with ",
                       toString(n+1), methodName,
                       " nodes", sep=""))
    points(z,f, col = 'red')
    legend('topright', legend=c('function','interpolation'),
           col = c('black','red'), lwd=1)
}
    
    if (Return){
        return(f)
}
}



piecewiseInterpol=function(n,nInt,a,b,neval, nodes = "equi", FUN, Plot){
## @param n : the degree of the interpolating polynomial on each
## subinterval
## @param nInt : the number of sub-intervals
## @param a, b : endpoints of the interval
## @param neval : the number of points on the interpolating grid (on
## each subinterval)
## @param nodes : string, either "equi" (default) for equidistant
## Lagrange interpolation (on each subinterval) or "cheby" for
## chebyshev nodes.
## @param FUN the function to be interpolated
## @param Plot : logical. Should the result be plotted ?
## @return : a matrix with 2 rows and neval * nInt -neval + 1:
## values of the interpolated funtion on a regular grid (first row)
## and the corresponding abscissas (second row).
    intEndPoints = seq(a,b,length.out = nInt+1)
    f = c()
    z = c()
    for (m in 1:nInt){
        A = intEndPoints[m]; B = intEndPoints[m+1]
        fm = interpolLagrange(n=n, a=A, b=B, neval=neval,  nodes = nodes, FUN=FUN, Plot=FALSE)
        zm = seq(A,B, length.out = neval)
        
        if( m >= 2){
        ## remove first element of zm, fm to avoid
        ## duplicate values of the interpolating vector
            zm = zm[-1]
            fm = fm[-1]
       
        }
        z = c(z,zm)
        f = c(f,fm)
    }

    if (Plot == 1){
        if (nodes == "equi") {methodName = " equidistant "}
        else {methodName = " Chebyshev "}

        plot(z, sapply(z,FUN),type="l")
        title(main = paste("Piecewise Lagrange interpolation with ",
                           toString(n+1), methodName, " nodes on ",
                           toString(nInt), " Intervals", sep=""))

        lines(z,f, col='red', lwd=2)
        
        legend('topright', legend = c('function','interpolation'),
               lwd=c(1,2), col=c('black','red'))

    }
    return(rbind(f,z) )
}

choose=function(neval, a, b, nodes = 'equi', FUN, Plot){
## Algorithme qui permet de choisir automatiquement la paire (M,n) 
## avec un budget total neval d'évaluations qui minimisent l'erreur d'interpolation.    
## @param neval : le budget maximal d'évaluations
## @param a : point d'extrémité gauche de l'intervalle total
## @param b : point d'extrémité droit de l'intervalle total
## @param nodes :string, soit "equi" (default) pour une interpolation 
## de Lagrange avec noeuds équidistants ou "Cheby" utilisant les noeuds 
## de Chebyshev.
## @param FUN: la fonction à interpoler 
## @param Plot : logical. Affiche (True) ou non (False) le graphique.
## @return : vecteur de taille 2 : (M,n) qui minimise l'erreur.

    normSupError = Inf
    z = c()
    error = c()
    n_int = 0
    
    # la différence entre "equi" et "cheby" est la façon de calculer neval_int
    
    if (nodes == "equi"){
        for (i in 1:(neval - 1)){
            n_int = n_int + 1 
            M_int = floor((neval - 1)/n_int)
            
            # Répartition d'environ 2000 points équirépartis sur [a,b]
            # neval_int est défini de telle manière car piecewise prend le nombre de points par sous intervalle
            neval_int = floor((2000-1)/M_int + 1)
            
            # la deuxième ligne de piecewiseInterpol renvoie les points utilisés pour l'interpolation.
            Fonction = FUN(piecewiseInterpol(n_int, M_int, a, b, neval_int, nodes, FUN, FALSE)[2,])
            Interpol = piecewiseInterpol(n_int, M_int, a, b, neval_int, nodes, FUN, FALSE)[1,]
            normSupError_int = max(abs(Fonction - Interpol))
                     
            z = c(z,i)
            error = c(error,normSupError_int)
                  
            if(normSupError_int<normSupError){
                n = n_int
                M = M_int
                normSupError = normSupError_int
            }
        }
    }
    else if (nodes == "cheby"){
        for (i in 1:(neval - 1)){
            n_int = n_int + 1
            M_int = floor(neval/(n_int+1))
            neval_int = floor((2000-1)/M_int + 1)
            
            Fonction = FUN(piecewiseInterpol(n_int, M_int, a, b, neval_int, nodes, FUN, FALSE)[2,])
            Interpol = piecewiseInterpol(n_int, M_int, a, b, neval_int, nodes, FUN, FALSE)[1,]
            normSupError_int = max(abs(Fonction - Interpol))
            z = c(z,i)
            error = c(error,normSupError_int)
                  
            if(normSupError_int<normSupError){
                n = n_int
                M = M_int
                normSupError = normSupError_int
            }
        }
    }
        
    else{stop("the nodes must be either 'equi' or 'cheby'") }


    if (Plot == 1){
        
        if (nodes == "equi") {methodName = " equidistant "}
        else {methodName = " Chebyshev "}

        plot(z, error,type="l")
        title(main = paste("Erreur en fonction du degre n avec ",
                           methodName, " nodes"))

    }
    
return (c(M,n, error[n],error[n+3]))  
}



trapezeInt =function(FUN,a,b,M){
    ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
    ##' @param FUN : the function to be integrated
    ##' @param a, b : interval end points 
    ##' @param M : number of intervals (each of size (b-a)/M)
    ##' @return: the value of the composite trapezoidal quadrature. 
    x = seq(a,b, length.out= M+1)
    y = sapply(x, FUN)
    h = (b-a)/M 
    q = h*(-0.5*y[1] + sum(y) - 0.5*y[length(y)])
    
    return(q)
}


refineTrapeze=function(FUN,a,b,M,q){
    ##' refinement of the subdivision step: incremental method
    ##' @param FUN : the function to be integrated
    ##' @param a, b : interval end points 
    ##' @param M : initial number of intervals (each of size (b-a)/M)
    ##'  having been used to compute q
    ##' @param  q : the value of the trapezoidal  quadrature method
    ##'  of stepsize (b-a)/M
    ##' @return : the value of the quadrature for a stepsize h' = h/2
    h = (b-a)/M
    x =  a + seq(1, 2*M-1, by=2)*h/2  
    ##  x : a vector of size M :
    ##     the additional abscissas where 'fun' must be evaluated.
    y = sapply(x, FUN)
    Q = 0.5*(q + h*sum(y))
    return(Q)
}



simpsonInt = function(FUN,a,b,M){
    ##' Simpson integration via trapeze rule
    ##' uses the fact that 
    ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
    h = (b-a)/M;
    qtrapeze = trapezeInt(FUN, a, b, M)  
    qrefined = refineTrapeze(FUN, a, b, M, qtrapeze)
    q =  4/3 * (qrefined - 1/4 * qtrapeze) 
    return(q)
}



 evalErrSimpson=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 ) 
  q = simps_h2  
  E = (simps_h - simps_h2)/15
    return(c(E,q))
}




richardson = function(FUN,n,t,delta){
    ## Calcule le tableau des differences  divisees en 0 du 
    ## polynome d'interpolation en t,delta t, ... delta^n t
    ## renvoie un vecteur de taille n+1:
    ## le vecteur des A_{k,k}, k= 0 .. n 
    ## (pas la matrice).   
    ## La meilleure approximation est le dernier element A[n+1].
    ##
    lx = log(t)  +  log(delta) *(0:n)
    x = exp(lx) 
    A = sapply(x,FUN) 
    for( j in 2:(n+1)){
        A[j : (n+1) ] =  (A[j:(n+1)] - delta^(j-1)* A[(j-1):n])/ (1 - delta^(j-1))
    }
    return(A)
}


romberg =function(FUN,n,a,b,M){## methode de Romberg avec n etapes
    ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
    ## pas initial h = (b-a)/M
    h= (b-a)/M 
    A = rep(0, n+1)
    A[1] = trapezeInt(FUN,a,b,M);
    Mc = M
    ## initialisation des differences divisees
    for( i in 2:(n+1)){
        A[i] = refineTrapeze( FUN,a,b, Mc, q= A[i-1])
        Mc = 2*Mc 
    }
    delta = 1/4;
    for (j in 2:(n+1)){
          A[j : (n+1) ] = (A[j:(n+1)] - delta^(j-1)* A[(j-1):n])/ (1 - delta^(j-1))
    }
    return(A)
}
