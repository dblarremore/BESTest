% Code for the Blockmodel Entropy Significance Test (BESTest)
% May 3, 2017
%
% Described in the paper:
% "The ground truth about metadata and community detection."
% Peel, Larremore, Clauset. Science Advances, 2017. 
% http://danlarremore.com/metadata 
% Comments or questions to larremore@santafe.edu
%
% Contents:
% 
% p = BESTest(adjMtx,partition,nSamples,modelName)
% 
% OUTPUT
%   * p - p value as described in the paper
% INPUTS
%   * adjMtx - NxN undirected network adjacency matrix. Code will check to ensure
%   that the matrix is either symmetric or triangular.
%   * partition - Nx1 or 1xN vector in which partition(i) is an integer
%   that enumerates which group vertex i belongs to. 
%   * nSamples - number of samples used to compute p. Recommended 10k or
%   higher for confident results. 
%   * modelName - string with desired model. Four options:
%       * 'SBMpoisson'
%       * 'dcSBMpoisson'
%       * 'SBMbernoulli'
%       * 'dcSBMmultinomial'

function p = BESTest(A,g,N,modelName)
% check inputs
isValidN(N);
isValidModelName(modelName);
A = isValidMtx(A);
g = isValidPartition(g);
h = zeros(N,1);
switch modelName
    case 'SBMpoisson'
        % note: multiplying by -1 for Likelihoods.
        [h_partition,~] = -1*SBMlikelihood(A,g);
        for n=1:N
            [h(n),~] = -1*SBMlikelihood(A,shuffle(g));
        end
    case 'dcSBMpoisson'
        % note: multiplying by -1 for Likelihoods.
        [~,h_partition] = -1*SBMlikelihood(A,g);
        for n=1:N
            [~,h(n)] = -1*SBMlikelihood(A,shuffle(g));
        end
    case 'SBMbernoulli'
        h_partition = SBMentropy(A,g);
        for n=1:N
            h(n) = SBMentropy(A,shuffle(g));
        end
    case 'dcSBMmultinomial'
        h_partition = DCSBMentropy(A,g);
        for n=1:N
            h(n) = DCSBMentropy(A,shuffle(g));
        end
end
p = sum(h <= h_partition)/n;
return
end

%%%%%%%%%% - Functions Called - %%%%%%%%%%

%%% Check to see if partition is valid
function [answer] = isValidPartition(g)
% integer?
if sum(mod(g,1)==0) ~= 0
    error('Partition vector is non-integer. Let g(i) = integer ID of vertex i`s group.')
end
% minimum group index
k = min(g);
% set minimum to 1
g = g - k + 1;
K1 = max(g);
K2 = length(unique(g));
if K1 ~= K2
    [ung,~,IC] = unique(g);
    for i=1:length(ung)
        h(find(IC==ung(i))) = i;
    end
    g = h;
end
return g
end

%%% Check to see if adjacency matrix is valid
function [answer] = isValidMtx(A)
if size(A,1) ~= size(A,2)
    error('Adjacency matrix not square');
else
isSymmetric = sum(sum(A~=A'))==0;
isUpperTriangular = sum(sum(A-triu(A)))==0;
isLowerTriangular = sum(sum(A-tril(A)))==0;
if sum(isSymmetric+isUpperTriangular+isLowerTriangular)==0
    error('Adjacency matrix is not symmetric or triangular')
end
if isSymmetric
    return
else
    A = A+A';
    return
end
end

%%% Check to see if N is valid
function [answer] = isValidN(N)
if N < 1 
    error('nSamples must be positive')
elseif mod(N,1)~=0
    error('nSamples must be integer')
end
answer = True;
return
end

%%% Check to see if modelName is valid
function [answer] = isValidModelName(modelName)
answer = True;
if ~ischar(modelName)
    error('modelName must be a character string')
end
if strcmp(modelName,'SBMpoisson')
    return
elseif strcmp(modelName,'dcSBMpoisson')
    return
elseif strcmp(modelName,'SBMbernoulli')
    return
elseif strcmp(modelName,'dcSBMmultinomial')
    return
else
    error('Invalid modelName. (SBMpoisson; dcSBMpoisson; SBMbernoulli; dcSBMmultinomial)')
end
end