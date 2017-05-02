# MATLAB for Blockmodel Entropy Significance Test (BESTest)
May 3, 2017
Described in the paper:
"The ground truth about metadata and community detection."
Peel, Larremore, Clauset. Science Advances, 2017. 
http://danlarremore.com/metadata 
Comments or questions to larremore@santafe.edu
# Usage:
BESTest.m is should be called 
p = BESTest(adjMtx,partition,nSamples,modelName)
# OUTPUT
  * p - p value as described in the paper
# INPUTS
  * adjMtx - NxN undirected network adjacency matrix. Code will check to ensure
  that the matrix is either symmetric or triangular.
  * partition - Nx1 or 1xN vector in which partition(i) is an integer
  that enumerates which group vertex i belongs to.
  * nSamples - number of samples used to compute p. Recommended 10k or
  higher for confident results.
  * modelName - string with desired model. Four options:
      * 'SBMpoisson'
      * 'dcSBMpoisson'
      * 'SBMbernoulli'
      * 'dcSBMmultinomial'
# Example Code:
see lazegaLawyersDemo.m for sample code and usage, fully reproducing Table 1 from the manuscript
