% Code for the Blockmodel Entropy Significance Test (BESTest)
% May 3, 2017
%
% Described in the paper:
% "The ground truth about metadata and community detection."
% Peel, Larremore, Clauset. Science Advances, 2017. 
% http://danlarremore.com/metadata
% 
% Comments or questions to larremore@santafe.edu

function [SBM,DCSBM] = SBMlikelihood(A,g)
% L(sbm,dcsbm) = sbmLikelihood(adjMtx,partition)
% This function evaluates the unnormalized log likelihood function L
% We calculate the Maximum Likelihood values of theta and omega, which can
% be calculated from the partition.

% number of nodes, N
N = size(A,1);
% degrees, k
k = sum(A);
% number of comms, K
K = max(g);
% degrees per group kappa
% nodes per group n
kappa = zeros(K,1);
n = zeros(K,1);
for i=1:K
    kappa(i) = sum(k(g==i));
    n(i)= sum(g==i);
end

% ML theta parameters (not actually used)
theta = zeros(N,1);
for i=1:N
    theta(i) = k(i)/kappa(g(i));
end

% calculate omega (equal to m)
[r,c,v] = find(A);
% let's be clever about how matlab constructs matrices
m = full(sparse(g(r),g(c),v,K,K));

intermediate = m.*log(diag(1./kappa)*m*diag(1./kappa));
intermediate(isnan(intermediate)) = 0;
DCSBM = sum(sum(intermediate));

intermediate = m.*log(diag(1./n)*m*diag(1./n));
intermediate(isnan(intermediate)) = 0;
SBM = sum(sum(intermediate));

return
end