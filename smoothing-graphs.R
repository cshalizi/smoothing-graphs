# Code for learning graphon models by smoothing adjacency matrices
# Cosma Shalizi, after conversation with Brian Karrer 8 May 2013
  # See end of file for license (2-claused simplified BSD license)


require(igraph)



# Learn a graphon model of a graph by smoothing and latent location
# Inputs: a graph adjacency matrix (net)
  # initial latent-location function for nodes (primer)
  # convergence scale (epsilon)
  # maximum iterations (maxit)
  # flag for iteratively adjusting the latent locations (iterative)
# Outputs: estimated function (graphon)
  # estimated latent locations (locs)
  # number of iterations used to convergence (iteration)
# TODO: Allow for priming with a vector and not just a function
# TODO: Abstract out to more general smoothers
  # Ideally, give the name of a smoother function as an argument, plus suitable
  # control settings (formula? bandwidth? what else?)
  # Make sure it returns a prediction function
graph.smoother <- function(net,primer,epsilon=0.01,maxit=100,iterative=TRUE) {
  require(mgcv) # for smoothing
  # Get an initial vector of locations
  v <- primer(net)
  # transform to the uniform(0,1) scale by ECDF
  u <- uniformize(v)
  # Set up the counters
  converged <- FALSE
  iteration <- 0
  # If we're not doing iterative approximation, reset no. of iterations
  if(!iterative) { maxit <- 1 }
  # until we converge or grow tired...
  while(!converged && (iteration < maxit)) {
    # Create a data-frame for smoothing
    z <- df.from.graph.and.latents(u,net)
    # Do the smoothing
      # spline smoothing on both variables at once, using logit link function
      # to ensure that estimated probabilities are in (0,1)
    graphon <- gam(a ~ s(u1,u2),data=z,family=binomial)
    # Get new latent locations, if we're into that
    if (iterative) {
       new.u <- update.u(net,u,graphon)
    } else {
       new.u <- u
    }
    # Test for convergence
    converged <- !iterative || (all(abs(new.u-u)<epsilon))
    # Update latent locations, iteration count
    u <- new.u
    iteration <- iteration+1
  }
  return(list(graphon=graphon,u=u,iteration=iteration))
}



# Prepare a data frame, suitable for smoothing, from a graph and a vector of
# latent node locations
  # The data frame should have the format of (latent1, latent2, adjacency)
# Inputs: vector of scalar node locations (u), graph on n nodes (graph)
# Output: a 3-column, n^2 row data frame
# Calls: adjacency.vector()
df.from.graph.and.latents <- function(u,net) {
  # Should have one latent location per node
  stopifnot(length(u) == vcount(net))
  # Need all pairwise combinations of latent locations
  u.pair <- expand.grid(u1=u,u2=u)
  return(data.frame(u.pair,a=adjacency.vector(net)))
}


# Improve latent location of nodes in a graph, given the graphon function, by
# maximum likelihood
# Inputs: network as an igraph object (net)
  # vector of latent locations (u)
  # graphon function (graphon)
# Outputs: new vector of latent locations
# Presumes: "graphon" really is a model object with a prediction method, and
  # looks for variables named "u1" and "u2"
update.u <- function(net,u,graphon) {
  # Check that u has the right length to go with net
  n <- length(u)
  stopifnot(n == vcount(net))
  # Create the adjacency vector, for likelihood purposes
  adj <- adjacency.vector(net)
  # Calculate negative log-likelihood as a function of u, fixing net and graphon
  negloglike <- function(u) { 
    # Get all the u-coordinate combinations
    u.pair <- expand.grid(u1=u,u2=u)
    # Get predicted edge probabilities
    p <- predict(graphon,newdata=u.pair,type="response")
    # Take either p or 1-p, depending on whether an edge exists
    l <- ifelse(adj>0,log(p),log(1-p))
    # Take logs, sum, take minus, return
    return(-sum(l))
  }
  # Optimize the log-likelihood via the Newton-style BFGS method
    # Presumes that the whole space is available rather than being [0,1]
  fit <- optim(par=u,fn=negloglike,method="L-BFGS-B",lower=rep(0,n),
     upper=rep(1,n),control=list(trace=1))
  # We just want the estimated locations
  return(fit$par)
}

# Sort an adjacency matrix by degree
# Inputs: network (net), represented by an igraph object
# Output: vector of ordered positions
degree.primer <- function(net) {
  return(order(degree(net)))
}

# Sort an adjacency matrix by spectral clustering
# Inputs: network as an igraph object (net)
# Outputs: vector of ordered positions
spectral.primer <- function(net) {
  laplacian <- graph.laplacian(net)
  # Need the next-to-zero eigenvector
    # TODO: Find a faster way to get the important eigenvector
  n <- vcount(net)
  spectral.vector <- eigen(laplacian)$vector[,n-1]
  return(spectral.vector)
}

# Turn a graph into an adjacency _vector_, i.e., an adjacency matrix in,
# to be definite, column-major order
# Inputs: a graph on n nodes, either in igraph format or a matrix (net)
# Outputs: a vector of length n^2
adjacency.vector <- function(net) {
  if(is.igraph(net)) {
    net <- get.adjacency(net)
  }
  return(as.vector(net))
}


# Transform arbitrary vector to uniform (0,1) scale via ECDF
  # HACK: the ECDF of the max is always 1.0, which is on the boundary;
  # multiply the ECDF values 0.999 to avoid this.
  # TODO: Figure out a more elegant approach
# Input: vector of numbers (x)
# Output: corresponding quantiles
uniformize <- function(x) {
  return(0.999*ecdf(x)(x))
}


# Transform arbitrary vector to standard Gaussian distribution via ECDF
  # HACK: the ECDF of the max is always 1.0, which is always Inf on a Gaussian
  # scale; multiply the ECDF by 0.999 to avoid this.
  # COMMENT: This didn't work very well in simulation trials, too many points
  # found ways to veer off to infinity to take advantage of ill-estimated
  # bits of the graphon  
# Input: vector of numbers (x)
# Output corresponding CDF-and-quantile-transformed vector
gaussianize <- function(x) {
  return(qnorm(uniformize(x)))
}



####################################
#      LICENSE                     #
####################################

# Copyright (c) 2013, Cosma Shalizi
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of their respective employers.
