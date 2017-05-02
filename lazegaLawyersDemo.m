% Code for the Blockmodel Entropy Significance Test (BESTest)
% May 3, 2017
%
% Described in the paper:
% "The ground truth about metadata and community detection."
% Peel, Larremore, Clauset. Science Advances, 2017. 
% http://danlarremore.com/metadata 
% Comments or questions to larremore@santafe.edu

% This code reproduces some results of the Blockmodel Entropy Sinificance 
% Test from the paper on the Lazega Lawyers dataset. 


% Load in Lazega Lawyers Data
% NB: adj_* are adjacency matrices. met_* are metadata vectors
load lazega.mat

% BESTest code is called as follows:
% p = BESTest(adjMtx,partition,nSamples,modelName)
nSamples = 10000;
modelName = 'SBMbernoulli';

fprintf('\ngetting p values for the friendship network')
adjMtx = adj_friend;
p(1,1) = BESTest(adjMtx,met_status,nSamples,modelName); fprintf('.');
p(1,2) = BESTest(adjMtx,met_gender,nSamples,modelName); fprintf('.');
p(1,3) = BESTest(adjMtx,met_office,nSamples,modelName); fprintf('.');
p(1,4) = BESTest(adjMtx,met_practice,nSamples,modelName); fprintf('.');
p(1,5) = BESTest(adjMtx,met_school,nSamples,modelName); fprintf('.\n');

fprintf('getting p values for the cowork network')
adjMtx = adj_work;
p(2,1) = BESTest(adjMtx,met_status,nSamples,modelName); fprintf('.');
p(2,2) = BESTest(adjMtx,met_gender,nSamples,modelName); fprintf('.');
p(2,3) = BESTest(adjMtx,met_office,nSamples,modelName); fprintf('.');
p(2,4) = BESTest(adjMtx,met_practice,nSamples,modelName); fprintf('.');
p(2,5) = BESTest(adjMtx,met_school,nSamples,modelName); fprintf('.\n');

fprintf('getting p values for the advice network')
adjMtx = adj_advice;
p(3,1) = BESTest(adjMtx,met_status,nSamples,modelName); fprintf('.');
p(3,2) = BESTest(adjMtx,met_gender,nSamples,modelName); fprintf('.');
p(3,3) = BESTest(adjMtx,met_office,nSamples,modelName); fprintf('.');
p(3,4) = BESTest(adjMtx,met_practice,nSamples,modelName); fprintf('.');
p(3,5) = BESTest(adjMtx,met_school,nSamples,modelName); fprintf('.\n\n');

% print results
fprintf('Net\tStatus\tGender\tOffice\tPract.\tSchool\n')
fprintf('Friend\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',p(1,1),p(1,2),p(1,3),p(1,4),p(1,5))
fprintf('Cowork\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',p(2,1),p(2,2),p(2,3),p(2,4),p(2,5))
fprintf('Advice\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',p(3,1),p(3,2),p(3,3),p(3,4),p(3,5))
