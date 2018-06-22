## About

This directory contains the scripts that were used to genereta and analyze the data presented in the manuscript. 
So they may serve as a tool to reproduce the data, and a reference for new studies.

The output of this scripts can also be examined in the respective $modcell2/problems folders.

## Computation reproducibility:

The MOEA algorithm relays on a random number generator to find solutions. 
When using the serial version, setting the state of the random number
generator before running the algorithm should be enough to obtain reproducible
results according to Matlab's documentation.
However, when using the parallel version, that method is no longer valid and currently there does not seem to be a general way to make parallel computations reproducible (Matlab has
an option for this but only for functions in the Statistics and Machine
Learning Toolbox). That said, the convergence criteria is robust enough so that the overall
solution (i.e. pareto front) looks the same in different runs.

When going through the examples, if you re-run the calculations (be aware that
this can take a while for the larger problems), you will obtain  similar
solutions but with different indices, and maybe different alternative soluitons, thus if some analysis script referes to
a specific index, they correspond to the particular output included in the output folder.
