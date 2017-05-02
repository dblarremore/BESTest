% Code for the Blockmodel Entropy Significance Test (BESTest)
% May 3, 2017
%
% Described in the paper:
% "The ground truth about metadata and community detection."
% Peel, Larremore, Clauset. Science Advances, 2017. 
% http://danlarremore.com/metadata
% 
% Comments or questions to larremore@santafe.edu

function [H,p] = DCSBMentropy(A,g)
% Adjacency matrix A
% metadata vector of group assignments g
K = max(g);
% # of vertices
N = size(A,1);

% count edges in and between each group, e_rs
for r=1:K
    for s=1:K
        e(r,s) = sum(sum(A(g==r,g==s)));
    end
end
m = sum(e);
for r=1:K
    for s=1:K
        e(r,s) = e(r,s)/(m(r)*m(s));
    end
end
% degrees, k
k = sum(A);
% pij = ki ers kj / er es 2m
% we've precomputed a bunch to speed this up. 
p = k'*k/sum(m);
% for i=1:N
%     for j=i:N
% %         p(i,j) = k(i)*k(j)*e(g(i),g(j)) / (sum(m)*m(g(i))*m(g(j)));
% %         p(j,i) = p(i,j);
%         p(i,j) = p(i,j) * e(g(i),g(j));
%     end
% end
p = p.*e(g,g);
H = multinomialEnt(triu(p),sum(m));
end

function h = multinomialEnt(x,n)
% h = entropy of a multinomial distribution with parameters x and n draws
% multinomial parameters, but only those that are > 0.
X = x(:);
X(X==0) = [];
% number of categories in the multinomial distribution, m
b = length(X);
%Eq (6), Jacek Cichon and Zbigniew Go?ebiewski, DMTCS Proc 2012. 
h = (1/2) * ( (b-1)*log(2*pi*n*exp(1)) + sum(log(X)) );
end