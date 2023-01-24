# Title: Knorr-Held Constraints in INLA
# Author: Taylor Okonek 

library(INLA)
library(tidyverse)

# If ohioadj.asc is not in your current working directory, change this to be the path 
# directing to ohioadj.asc
ohio_graph_path <- "ohioadj.asc"


# Data set-up -------------------------------------------------------------

# read data frame
data.inla <- read.table("http://faculty.washington.edu/jonno/SISMIDmaterial/ohiodata2_spacetime.txt",header=T)
head(data.inla)

#--- add all the necessary columns to the data ---#
# Number of areas and number of time points, for later use
N <- 88 # number of areas
S <- 21 # number of time points
ohiosp <- data.frame(fips=data.inla$county,
                     county= data.inla$county, 
                     year=data.inla$time,
                     expected=data.inla$E,
                     deaths=data.inla$Y)

# create copies of the time indicator for a structured and unstructured random effect, if needed
ohiosp$time.unstruct <- ohiosp$time.struct <- ohiosp$year

# read in and create ICAR adjacency matrix
neighbors <- read.table(ohio_graph_path, header=T)
m <- neighbors$num
# alternatively, could read in a .graph file, so if you have a .graph file this set-up
# will probably look a bit different
ICAR_prec <- diag(N)
diag(ICAR_prec) <- m
#
for(i in 1:N){ 
  neighs <- neighbors[i,3:ncol(neighbors)][neighbors[i,3:ncol(neighbors)]!=0]
  ICAR_prec[i,neighs] <- (-1)
}

# quick check: If you coded your ICAR precision matrix correctly, both rows and
# columns will sum to zero
apply(ICAR_prec,1,sum)
apply(ICAR_prec,2,sum)


# Main effects model ------------------------------------------------------

# inla.rw is a function that will create the precision matrix of a RW(1 or 2)
# so that we can use the bym2 parameterization for the temporal component

# we'll use a RW1 in this example
# Note: if you have a RW2 you *may* need additional constraints, depending on what else is in your model
RW_prec <- INLA:::inla.rw(n = S, order = 1, scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                   sparse = TRUE)

# define a baseline formula with no interaction terms
# + BYM2 in space
# + BYM2 in time (RW1)
formula_base <- deaths ~ 1 + 
  f(year, model = "bym2", 
    scale.model = TRUE,
    constr = TRUE,
    rankdef = 1, # good practice to specify rank deficiency of precision matrix
    graph = RW_prec,
    hyper = list(
      phi = list(
        prior = "pc",
        param = c(0.5, 2/3)),
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01)))) +
  f(county, model = "bym2", 
    scale.model = TRUE,
    constr = TRUE,
    rankdef = 1, # good practice to specify rank deficiency of precision matrix
    graph = ICAR_prec,
    hyper = list(
      phi = list(
        prior = "pc",
        param = c(0.5, 2/3)),
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01)))) 

# Takes a bit of time to run
mod_main <- inla(formula_base,
             data = ohiosp,
             family = "poisson",
             E = expected,
             control.predictor = list(compute = TRUE),
             control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))

# Check the number of constraints in our model
# We should have:
# - 1 sum-to-zero constraint for the BYM2 in space
# - 1 sum-to-zero constraint for the BYM2 in time
# - 2 constraints total

mod_main$misc$configs$constr$nc

# Type I Interaction ------------------------------------------------------

# Type I Interaction: 
# - IID Space x IID Time 
# - No constraints

# create an ID column for the space-time interaction
# a unique number for each time.space combination (88*21 = 1848)
ohiosp$spacetime.unstruct<- 1:(N *S)


# NOTE: you need to be careful when setting the rankdef parameter if your graph has
# more than one connected component! For the temporal term this shouldn't be a problem,
# but if your graph has islands you need to be careful

# Add Type I Interaction term to base formula
formulaI <- update(formula_base, ~. + f(spacetime.unstruct, model="iid",
                                         hyper = list(prec = list(param = c(1, 0.01)))))

modI <- inla(formulaI,
             data = ohiosp,
             family = "poisson",
             E = expected)

# Check the number of constraints in our model
# We should have:
# - 1 sum-to-zero constraint for the BYM2 in space
# - 1 sum-to-zero constraint for the BYM2 in time
# - 2 constraints total

modI$misc$configs$constr$nc

# Type II Interaction ----------------------------------------------------------

# Type II Interaction: 
# - IID Space x RW(1) Time 
# - Constraints: The RW1 in each area needs to sum to 1
# - Number of constraints = rank deficiency of precision = N

# creating constraint matrix needed for Type II Interaction
A <- matrix(0, nrow = N, ncol = N * S)
for (i in 1:N) { 
  A[i, which((1:(N * S))%%N == i - 1)] <- 1
}

# specify constraints in INLA-ready format
constr.st <- list(A = A, e = rep(0, dim(A)[1]))

### defining Kronecker product for the Type II Interaction

# Kronecker product between RW1 and IID space term
# order matters here! In our data set, the time ordering take precedent over the space ordering
# so we must have the RW on the left and space on the right
# Optionally, we scale the matrix beforehand as well
scaled_RW_prec <- inla.scale.model(RW_prec, list(A = matrix(1, 1, dim(RW_prec)[1]), e = 0))
R <- scaled_RW_prec %x% diag(N)

formulaII <- update(formula_base, ~. + f(spacetime.unstruct, 
                                      model = "generic0", Cmatrix = R, extraconstr = constr.st, 
                                      rankdef = N, 
                                      param = c(1, 0.01)))

### Time difference of 15.01259 mins on my machine
ptm <- Sys.time()
# modII <- inla(formulaII,
#              data = ohiosp,
#              family = "poisson",
#              E = expected)
Sys.time() - ptm

# Check the number of constraints in our model
# We should have:
# - 1 sum-to-zero constraint for the BYM2 in space
# - 1 sum-to-zero constraint for the BYM2 in time
# - N constraints for the Type II Interaction
# - 2 + N constraints total

modII$misc$configs$constr$nc

# Type III Interaction ----------------------------------------------------------

# Type II Interaction: 
# - ICAR Space x IID Time 
# - Constraints: The ICAR at each time needs to sum to 1
# - Number of constraints = rank deficiency of precision = S

### creating constraint needed for Type III Interaction
A <- matrix(0, nrow = S, ncol = N * S)
for (i in 1:S) {
  # The ICAR at each time point needs to sum to 1
  A[i, ((i - 1) * N + 1):(i * N)] <- 1
}

### defining Kronecker product for the Type III Interaction
# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_prec, 
                                           constr = list(A = matrix(1, 1, dim(ICAR_prec)[1]), e = 0))

# Kronecker product between IID x ICAR
R <- diag(S) %x% ICAR_prec 

# specify constraints in INLA-ready format
constr.st <- list(A = A, e = rep(0, dim(A)[1]))

formulaIII <- update(formula_base, ~. + f(spacetime.unstruct, 
                                      model = "generic0", Cmatrix = R, extraconstr = constr.st, 
                                      rankdef = S, 
                                      param = c(1, 0.01)))

### Time difference of 1.319204 mins on my personal machine
ptm <- Sys.time()
# modIII <- inla(formulaIII,
#               data = ohiosp,
#               family = "poisson",
#               E = expected)
Sys.time() - ptm


# Check the number of constraints in our model
# We should have:
# - 1 sum-to-zero constraint for the BYM2 in space
# - 1 sum-to-zero constraint for the BYM2 in time
# - S constraints for the Type III Interaction
# - 2 + S constraints total

modIII$misc$configs$constr$nc


# Type IV Interaction ----------------------------------------------------------

# Specifying constraints needed for Type IV Interaction
time_constr <- matrix(0, N, N * S)
for (i in 1:N) {
  time_constr[i, which((1:(N * S))%%N == i - 1)] <- 1
}
space_constr <- matrix(0, S-1, N * S)
# theoretically it's enough to only go through N-1 here
# note that we could have instead done 1:(S-1) in the first loop and 1:N in the second
for (i in 1:(S-1)) { 
  space_constr[i, ((i - 1) * N + 1):(i * N)] <- 1
}

tmp <- rbind(time_constr, space_constr) 
constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))

# Kronecker product between ICAR x RW1
R <- scaled_RW_prec %x% scaled_ICAR_prec

formulaIV <- update(formula_base, ~. + f(spacetime.unstruct, 
                                  model = "generic0", Cmatrix = R, extraconstr = constr.st,
                                  rankdef = N + S - 1, 
                                  param = c(1, 0.01)))

#  Time difference of ~ 20 minutes on my personal machine
ptm <- Sys.time()
# modIV <- inla(formulaIV,
#                data = ohiosp,
#                family = "poisson",
#                E = expected)
Sys.time()-ptm




# Check the number of constraints in our model
# We should have:
# - 1 sum-to-zero constraint for the BYM2 in space
# - 1 sum-to-zero constraint for the BYM2 in time
# - N + S - 1 constraints for the Type IV Interaction
# - N + S + 1 constraints total
modIV$misc$configs$constr$nc

# Constraints for a generic GMRF  -----------------------------------------------

# All of the above examples are for specific Type II - Type IV Knorr-Held models,
# assuming a RW1 for time. If you have some other GMRF that needs constraining (or you
# just want to go about this all in a different way that's perhaps intuitive if you like
# linear algebra)... you can do the following.

# If you want to FULLY constrain a matrix (so that is invertible and therefore ALL null 
# space is "removed"), you can construct the constraint matrix using eigenvectors, to enforce
# the constraint(s) Ax = e (where typically e = 0).

# Specificially, the null space of a matrix is spanned by the eigenvectors that correspond to
# the zero eigenvalues of a matrix. If a matrix has no zero eigenvalues, it is not improper
# and, therefore, invertible. This is why constraints work: we essentially "remove" the null
# space of a matrix by constraining it away.

# Let's take the RW1 as an example. 
# We can obtain the eigenvalues and eigenvectors of the matrix:

RW_eigen <- eigen(RW_prec)
RW_eigen$values

# As you can see above, there is one eigenvalue == 0 for this matrix. We can look at the 
# corresponding eigenvector:

RW_eigen$vectors[,21]

# The eigenvector(s) that correspond to the zero eigenvalues will make up the row(s) of our
# constraint matrix. For the RW, we do the following, since there's only one constraint:

A <- matrix(RW_eigen$vectors[,21], nrow = 1)

# Fun fact (that should hopefully make sense): We know that when we have an intercept in our
# model as well as a RW1, we need to put a sum-to-zero constraint on the RW1. Note that
# Ax in this case is sum(c * x), since it's just a one-row matrix of the same constant repeated.
# ...and sum(c * x = 0) is the same constraint as sum(x) = 0. So it makes sense!


# As an aside, an easy way to check if all of the constraints you created are actually spanning the
# null space of the precision matrix is just to multiply Ax and see if you get e (again,
# typically zero), as intended

A %*% RW_prec




