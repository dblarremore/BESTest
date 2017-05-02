% Code for the Blockmodel Entropy Significance Test (BESTest)
% May 3, 2017
%
% Described in the paper:
% "The ground truth about metadata and community detection."
% Peel, Larremore, Clauset. Science Advances, 2017. 
% http://danlarremore.com/metadata
% 
% Comments or questions to larremore@santafe.edu

function H = SBMentropy(A,pi)
% H = SBMentropy(A,pi)
% A is the adjacency matrix
% pi is a vector of integer group memberships

% Number of vertices N
N = size(A,1);
% Get the list of unique group numbers, q
q = unique(pi);
% Count the number of members of each group, n
for i=1:length(q)
    n(i) = sum(pi==q(i));
end

%%% NOTE ABOUT SELF LINKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the following, we assume that there are self-edges possible in the
% network. Moreover, we take the convention that IF there is an edge from 
% node i to itself, then A(i,i) = 2. This is inspected, then enforced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

self = diag(A);
if sum(self)==0
    hasSelfLoops = 0;
else
    if sum(self==1) > 0 && sum(self==2) > 0
        error('Self-loops have values of 1 AND 2 in your adjacency matrix.');
    elseif sum(self==1) > 0 && sum(self==2) == 0
        % We have self loops counted as 1 edge... modifying.
        hasSelfLoops = 1;
        A(1:end+1:end) = 2*A(1:end+1:end);
    elseif sum(self==1) == 0 && sum(self==2) > 0
        hasSelfLoops = 1;
    end
end

% Compute the group connectivity matrix, w
% Compute the group entropies, ent
w = zeros(length(q));
ent = zeros(length(q));
for r=1:length(q)
    for s=1:length(q)
        if r~=s
            w(r,s) = sum(sum(A(pi==q(r),pi==q(s)))) / (n(r)*n(s));     
            ent(r,s) = n(r)*n(s)*bernEnt(w(r,s));
        else
            if hasSelfLoops==1
                w(r,r) = sum(sum(A(pi==q(r),pi==q(r)))) / (n(r)*n(r));
                ent(r,r) = n(r)*n(r)*bernEnt(w(r,r));
            else
                w(r,r) = sum(sum(A(pi==q(r),pi==q(r)))) / (n(r)*n(r)-n(r));
                ent(r,r) = n(r)*(n(r)-1)*bernEnt(w(r,r));
            end
        end
    end
end

% Compute the final entropy of the ensemble, H
H = (1/2) * sum(sum(ent));
end

function h = bernEnt(x)
% h = entropy of a bernoulli process with parameter x
% This function uses log2, not log. Answer is in bits.
if x==1 || x==0
    % Force entropy = 0 when x={0,1}
    h=0;
    return
else
    h = -x * log2(x) - (1-x) * log2(1-x);
    return
end
end