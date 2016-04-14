function [ Hubs ] = hubCaps( Metrics )
%HUBCAPS Calculates group or individual hub measures
%
%   [Hubs] = hubCaps(Metrics);
% 
%   Binary - also weighted version
%
%   Total of 7 nodal measures
%   Runs best with Metrics = myMeasures(CIJ);
%
%   Inputs:     Metrics,    structure with values 
%
%   Outputs:    Hubs,       structure of hubs (overall, individual
%                           subject, metric specific)
%
% Michael Hart, University of Cambridge, December 2015

%% Initialise

%define inputs
deg = Metrics.degree; %degree centrality
cl = Metrics.closeness; %closeness centrality
bc = Metrics.node_betweeness; %node betweenness centrality
z = Metrics.zscore; %intra-module centrality
p = Metrics.participation; %inter-module centrality
v = Metrics.eigenvector; %eigenvector centrality
pr = Metrics.pagerank; %pagerank centrality

%basic parameters
m = size(deg,1); %number of nodes
n = size(deg,2); %1 if group, >1 if individuals: use for loop below

%check if binary
if max(rem(deg,1)) == 0
    disp('input is binary')
else
    error('input is weighted - please use weighted version')
end

%check if individual matrix or group array
if n == 1
    disp('input is single group average')
else
    disp('input is an array of individuals')
end

%initialise outputs
hub_rank = zeros(m,1); %vector of node with overall hub ranking
hub_score = zeros(m,n); %matrix of individual hubs

%individual metrics vectors
degHub = zeros(m,1); clHub = zeros(m,1); bcHub = zeros(m,1); 
zHub = zeros(m,1); pHub = zeros(m,1); vHub = zeros(m,1); 
prHub = zeros(m,1);

%% Calculate measures

for ii = 1:n; %per subject
    [~, I] = sort(deg(:,ii), 'descend'); %degree centrality
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1; %overall
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1; %individual
    degHub(I(1:10)) = degHub(I(1:10)) + 1;%metric

    [~, I] = sort(cl(:,ii), 'descend'); %closeness
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    prHub(I(1:10)) = clHub(I(1:10)) + 1;
    
    [~, I] = sort(bc(:,ii), 'descend'); %node betweenness centrality
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    bcHub(I(1:10)) = bcHub(I(1:10)) + 1;

    [~, I] = sort(z(:,ii), 'descend'); %z-score
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    zHub(I(1:10)) = zHub(I(1:10)) + 1;
    
    [~, I] = sort(p(:,ii), 'descend'); %participation co-efficient
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    pHub(I(1:10)) = pHub(I(1:10)) + 1;
    
    [~, I] = sort(v(:,ii), 'descend'); %eigenvector
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    vHub(I(1:10)) = vHub(I(1:10)) + 1;
    
    [~, I] = sort(pr(:,ii), 'descend'); %pagerank
    hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
    hub_score(I(1:10),ii) = hub_score(I(1:10),ii) + 1;
    prHub(I(1:10)) = prHub(I(1:10)) + 1;

end 

%% Parse outputs

Hubs.overall = hub_rank;
Hubs.individual = hub_score;
Hubs.metrics = [degHub clHub bcHub zHub pHub vHub prHub];

end

