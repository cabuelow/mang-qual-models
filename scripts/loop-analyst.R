library(LoopAnalyst)
library(QPress)
source('scripts/models.R')
source('scripts/helpers.R')

# get adjacency matrix

matrix <- adjacency.matrix(modelB)

# get the adjoint matrix

system.time(matrix.adj <- make.adjoint(matrix))

system.time(make.wfm(matrix.adj))

data(submatrix)
enumerate.loops(submatrix) # get all of the loops
feedback(submatrix)

data(cm.levins)
enumerate.paths(cm.levins,2,4)

## compute adjoint of a negative of a community matrix
# the value of a given element indicates the sum of feedback at all levels,
# under the assumption that, in the absence of relative specificiation
# of community matrix linkages, each level of feedback 
data(cm.dambacher)
make.adjoint(cm.dambacher)

feedback(cm.dambacher)

# absolute feedback

make.T(cm.dambacher)

## compute weighted feedback matrix

make.wfm(cm.dambacher)

