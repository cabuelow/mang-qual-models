library(LoopAnalyst)
library(QPress)

data(submatrix)
enumerate.loops(submatrix) # get all of the loops
feedback(submatrix)

data(cm.levins)
enumerate.paths(cm.levins,2,4)

## compute adjoint of a negative of a community matrix
data(cm.dambacher)
make.adjoint(cm.dambacher)

## compute weighted feedback matrix
data(cm.dambacher)
make.wfm(cm.dambacher)

