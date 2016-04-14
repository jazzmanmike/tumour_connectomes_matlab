function [ HeavyHubs ] = hubCapsHeavyTwo( normMetrics  )
%HUBCAPSHEAVYTWO Calculates group or individual hub measures
%
%   HeavyHubs = hubCapsHeavyTwo(normMetrics);
%
%   Inputs:     normMetrics, matrix with values 
%
%   Outputs:    Hubs,       structure of hubs (overall, individual
%                           subject, metric specific)
%
% Michael Hart, University of Cambridge, February 2016

%% Initialise

%define inputs
st = normMetrics(:, 1); %strength
bc = normMetrics(:, 2); %node betweenness centrality
z = normMetrics(:, 3); %intra-module centrality
p = normMetrics(:, 4); %inter-module centrality
v = normMetrics(:, 5); %eigenvector centrality

%basic parameters
m = size(normMetrics,1); %number of nodes

%initialise outputs
hub_rank = zeros(m,1); %vector of node with overall hub ranking
%individual metrics vectors
stHub = zeros(m,1); bcHub = zeros(m,1); zHub = zeros(m,1); 
pHub = zeros(m,1); vHub = zeros(m,1); 

%% Calculate measures

[~, I] = sort(st(:), 'descend'); %strength
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
stHub(I(1:10)) = stHub(I(1:10)) + 1;

[~, I] = sort(bc(:), 'descend'); %node betweenness centrality
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
bcHub(I(1:10)) = bcHub(I(1:10)) + 1;

[~, I] = sort(z(:), 'descend'); %z-score
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
zHub(I(1:10)) = zHub(I(1:10)) + 1;

[~, I] = sort(p(:), 'descend'); %participation co-efficient
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
pHub(I(1:10)) = pHub(I(1:10)) + 1;

[~, I] = sort(v(:), 'descend'); %eigenvector
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
vHub(I(1:10)) = vHub(I(1:10)) + 1;

%% Parse outputs

HeavyHubs.overall = hub_rank;
HeavyHubs.metrics = [stHub bcHub zHub pHub vHub];

end

