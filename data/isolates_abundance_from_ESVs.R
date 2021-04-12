### Juan Diaz-Colunga
### Apr 11 2021

### ----------------------------------------------------------------------
### This script shows an example and defines a function to obtain the
### species composition of a sample from the sequencing data when there is
### noise introduced by e.g. sequencing error.
### ----------------------------------------------------------------------

library(matlib)
library(nloptr)

### This method has limitations. First, we have no way to know how many copies of the
### 16S rRNA gene each species has, and we know this copy number is variable across
### species. The only way to account for this would be through full genome sequencing
### of all the species we are interested in (ideally every species from every community,
### which is not viable in our setup). We are forced to make the assumption that all
### species have the same, or at least very similar, number of copies of the 16S. If
### this is not the case, our method could be innacurate; however this limitation is
### inherent to the 16S sequencing technique and falls out of our scope for this work.
### The second limitation is the fact that we do not know how every species in our
### communities maps to the ESV space (even in relative terms), i.e., the fraction of
### sequences of each type that are found when each individual species is sequenced
### separately. We do have this information for the dominants of each community (which
### were isolated and sequenced individually), and in practice the abundances of these
### dominants are what we are interested in. We can bypass this limitation by assuming
### that every ESV that we find in a community that does not match one of our dominants
### corresponds unequivocally to a single species. This need not be the case in the most
### general scenario (many species can share the same ESV, and a single species may have
### several copies of the 16S each matching different ESVs). However, this assumption
### will only hamper accuracy in the abundance of sub-dominant species, which we are not
### interested in for this specific application.
### 
### Q: matrix mapping isolates to sequenced ESVs. Each column corresponds to a species,
### each row corresponds to an ESV. The element (i,j) of the matrix represents the
### fraction of sequences that match the i-th ESV when species j is sequenced.

# mock species-to-ESVs matrix
Q <- matrix(c(1,0,0,0,0,
              0,0.75,0,0.25,0,
              0,0.9,0,0,0.1,
              0,0,1,0,0),
            nrow=5,
            byrow=F)
Q <- rbind(Q,matrix(0,ncol=ncol(Q),nrow=5))
colnames(Q) <- paste('isolate_',1:ncol(Q),sep='')
rownames(Q) <- paste('seq_',1:nrow(Q))

# add as many pseudo-species as necessary to cover all ESVs (see note above)
n <- rownames(Q)[rowSums(Q) == 0]
Q <- cbind(Q,matrix(0,nrow=nrow(Q),ncol=length(n)))
for (i in 1:length(n)) {
  colnames(Q)[ncol(Q) - length(n) + i] <- paste('pseu_',i,sep='')
  Q[n[i],ncol(Q) - length(n) + i] <- 1
}

# mock community composition
set.seed(0)
cc <- matrix(runif(ncol(Q)),ncol=1)
rownames(cc) <- colnames(Q)
cc <- cc/sum(cc)

# mock sequencing data
seq <- Q %*% cc

# add noise to sequencing data
seq_noise <- seq + (0.3)*matrix(runif(nrow(Q)),ncol=1)
seq_noise[seq_noise<0] <- 0
seq_noise <- seq_noise/sum(seq_noise)

### If there was no noise (sequencing error etc.) then it would be straightforward to
### obtain community composition from the sequencing data
cc_inferred_1 <- Ginv(Q) %*% seq
rownames(cc_inferred_1) <- colnames(Q)

### But if there is noise this method can lead to problematic situations (e.g. negative
### species abundances)
cc_inferred_2 <- Ginv(Q) %*% seq_noise
rownames(cc_inferred_2) <- colnames(Q)

### To address this, we consider that our deconvolution process (going from sequence
### abundances to species abundances) gives us an array of the true community composition
### plus an array of noise:
###    Ginv(Q) %*% seq_noise = cc_true + error
### Where cc_true staisfies all(cc_true)>0 and sum(cc_true)==1.
### We want to obtain the cc_true vector that stasfies those conditions while simultaneously
### minimizing the error, namely sum(error^2). For that we will use the nloptr package.
### Thus, we want to minimize:
###    sum((Ginv(Q) %*% seq_noise - cc_true)^2)
### While imposing:
###    0 <= cc_true[1] <= 1
###    0 <= cc_true[2] <= 1
###    0 <= cc_true[3] <= 1
###    ...
### That is:
###    all(cc_true) >= 0
###    all(cc_true) <= 1
###
### We start by defining the function to minimize as a function of cc_true (here denoted as x)
### and its gradient
eval_f <- function(x) sum((Ginv(Q) %*% seq_noise - x)^2)
eval_grad_f <- function(x) -2*(Ginv(Q) %*% seq_noise - x)

### Our initial 'guess' would just be the result of doing Ginv(Q) %*% seq_noise, correcting
### negative values and mantaining normalization
x0 <- Ginv(Q) %*% seq_noise
x0[x0<0] <- 0
x0 <- x0/sum(x0)

### Constraint function (in order for x to be normalized, sum(x) must be 1 and thus this function
### must be 0)
eval_h <- function(x) sum(x) - 1
eval_jac_h <- function(x) rep(1,length(x))

### Options to pass to nloptr:
local_opts <- list('algorithm' = 'NLOPT_LD_MMA',
                   'xtol_rel' = 1.0e-7)
opts <- list('algorithm' = 'NLOPT_LD_AUGLAG',
             'xtol_rel' = 1.0e-8,
             'maxeval' = 1000,
             local_opts = local_opts,
             'check_derivatives' = FALSE,
             'check_derivatives_print' = 'all',
             'print_level' = 0)

### Invoke nloptr (lb and ub are the lower and upper bounds respectively, in this case 0 and 1
### for all the components of our vector x)
res0 <- nloptr(x0=x0,
               eval_f = eval_f,
               eval_grad_f = eval_grad_f,
               lb = rep(0,ncol(Q)),
               ub = rep(1,ncol(Q)),
               eval_g_eq = eval_h,
               eval_jac_g_eq = eval_jac_h,
               opts = opts)

### Now we have a solution for the community composition that is normalized, satisfies all species
### having non-negative abundances, and minimizes the residues:
cc_inferred_3 <- matrix(res0$solution,ncol=1)
rownames(cc_inferred_3) <- colnames(Q)

### At this point we have:
###    cc_inferred_1: the true community composition, resulting from deconvolutiong when there is no noise
###    cc_inferred_2: the unconstrained community composition resulting from a raw deconvolution with noise
###    cc_inferred_3: the constrained community composition resulting from deconvoluting while allowing for
###                   noise correction, mantaining normalization and non-negative abundances

### We can put all this in the form of a function to which we can pass the matrix Q and the sequencing data
### to give us the inferred community composition in return

# remove all variables from the example above
rm(list=ls())

# define function
species_composition_from_sequencing <- function(Q,seq,remove_unused_seqs=TRUE) {
  
  # inputs as matrices
  Q <- as.matrix(Q)
  seq <- as.matrix(seq)
  
  # remove sequences that are not represented
  if (remove_unused_seqs) {
    n <- rowSums(Q) >0 | seq>0
    Q <- Q[n,]
    seq <- seq[n,]
  }
  
  # add pseudo-species corresponding to sequences not in Q
  n <- rownames(Q)[rowSums(Q) == 0]
  if (length(n)) {
    Q <- cbind(Q,matrix(0,nrow=nrow(Q),ncol=length(n)))
    for (i in 1:length(n)) {
      colnames(Q)[ncol(Q) - length(n) + i] <- paste('pseu_',i,sep='')
      Q[n[i],ncol(Q) - length(n) + i] <- 1
    }
  }
  
  # define function to minimize and constraint function
  eval_f <- function(x) sum((Ginv(Q) %*% seq - x)^2)
  eval_grad_f <- function(x) -2*(Ginv(Q) %*% seq - x)
  eval_h <- function(x) sum(x) - 1
  eval_jac_h <- function(x) rep(1,length(x))
  
  # initial guess
  x0 <- Ginv(Q) %*% seq
  x0[x0<0] <- 0
  x0 <- x0/sum(x0)
  
  # options to pass to nloptr
  local_opts <- list('algorithm' = 'NLOPT_LD_MMA',
                     'xtol_rel' = 1.0e-7)
  opts <- list('algorithm' = 'NLOPT_LD_AUGLAG',
               'xtol_rel' = 1.0e-8,
               'maxeval' = 1000,
               local_opts = local_opts,
               'check_derivatives' = FALSE,
               'check_derivatives_print' = 'all',
               'print_level' = 0)
  
  # solve optimization problem
  res <- nloptr(x0=x0,
                eval_f = eval_f,
                eval_grad_f = eval_grad_f,
                lb = rep(0,ncol(Q)),
                ub = rep(1,ncol(Q)),
                eval_g_eq = eval_h,
                eval_jac_g_eq = eval_jac_h,
                opts = opts)
  
  # return solution
  species_composition <- matrix(res$solution,ncol=1,)
  rownames(species_composition) <- colnames(Q)
  
  return(species_composition)
  
}





