## Transition Matrix approach to get steady state distribution for Markov Chain

## Create an initial transition matrix. This is relatively arbitrary,
##  so you could play around by creating other initial matrices, and see
##  how different the final matrices look! (just make sure that the
##  rows of your matrix sum to 1)
P.mat = matrix(c(.1, .1, .3, .3, .2, .2, .4, .4, .3, .3, .2, .2, .4, .4, .1, .1), ncol = 4)

## Now specify your desired steady state probabilities. These are completely
##  up to you, so long as there is one for each state, and they sum to 1
pi.vec = c(.2, .3, .4, .1)

## Now I need to create the acceptance probability matrix filled with alpha
##  values. There are probably efficient ways to do this, but I think I will
##  use loops, because it is generally easier to see what is going on with loops.

# Create an empty matrix with the same dimensions as the P matrix
alpha.mat = matrix(NA, nrow = 4, ncol = 4)

# Now fill the off-diagonal elements of the alpha matrix
for (i in 1:4) {
  for (j in 1:4) {
    if (i != j) {  # remember, this is for off-diagonal probabilities only ...
      alpha.mat[i, j] = min((pi.vec[j] * P.mat[j, i]) / (pi.vec[i] * P.mat[i, j]), 1)
    }
  }
}

## Multiply each element in the P matrix by each element in the alpha matrix.
##  Note: it doesn't matter that the diagonal elements become NA - we will
##  replace these soon.
Pprime.mat = alpha.mat * P.mat   # Note: this is element multiplication, NOT matrix multiplication

## Now fill in the diagonals of P' by ensuring each row sums to 1
for (i in 1:4) {
  Pprime.mat[i, i] = 1 - sum(Pprime.mat[i, ], na.rm = TRUE)
}
Pprime.mat

## We can check this matrix's steady state distribution by taking the
##  first eigenvector of the transpose of this matrix, and scaling its
##  elements to sum to 1
my.pi.vec = eigen(t(Pprime.mat))$vector[, 1]
my.pi.vec = my.pi.vec / sum(my.pi.vec)

## and check
my.pi.vec
pi.vec

## Another way to check the steady state of the matrix is to take the
##  matrix to a large power i.e. assess the occupation probabilities
##  after a large number of steps. To do this, you will need to have an
##  easy way of taking matrices to a particular power - use the
##  matpower.func() function I have written below, unless you know of a
##  better way.

#### matpower.func ###########################################
## A function for taking a square matrix mat to a power pow ##
## Written by Steven Miller, some time in the distant past  ##
##############################################################
matpower.func = function (mat, pow) {
    if (nrow(mat) != ncol(mat)) stop("Matrix mat is not square!")
    symmetric = all(mat == t(mat))
    eig = eigen(mat)
    if (symmetric) return(eig$vec %*% diag(eig$val^pow) %*% t(eig$vec))
    else return(eig$vec %*% diag(eig$val^pow) %*% solve(eig$vec))
}

## For the modified transition matrix created from the original matrix
##  in the slides, the steady state is reached (at least, to 8 decimal places)
##  after 15 steps. If you started with your own transition matrix, you
##  might need more or fewer steps to reach steady state. I chosen 100 steps
##  here to be safe, although you might require even more!
matpower.func(Pprime.mat, 100)
pi.vec


## And finally, if you don't trust the mathematical theory, you could always
##  simulate 10,000 runs, each of 100 steps, through this Markov Chain, and
##  record the proportion of times the runs end up in each state.
result.vec = c()
for (i in 1:10000) {
  state = sample(1:4, 1)  # Randomly choose the starting state - theory says this shouldn't matter anyway, at steady state
  for (n in 1:100) {  # You would need to increase this if it takes longer to reach steady state for your transition matrix
    state = sample(1:4, 1, prob = Pprime.mat[state, ])  # Choose a new state at random, based on the row of transition probabilities for the state you are currently in
  }
  result.vec[i] = state
}
table(result.vec) / 10000
pi.vec

