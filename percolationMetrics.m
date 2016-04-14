function [ nodalMetrics, globalMetrics, hubScore ] = percolationMetrics( CIJ, n )
%PERCOLATIONMETRICS Summary of this function goes here
%   Used in nodeAttack, randomError, and focusedAttack
%   Adaption of myHeavyMeasures & hubCapsHeavy plus others
%   Therefore, designed to be 'silent'
%
%   [nodalMetrics, globalMetrics, hubScore] = percolationMetrics(CIJ, n);
%
%   Inputs:     CIJ,                single subject weighted matrix
%               n,                  number of removed nodes (optional)
%
%   Outputs:    nodalMetrics,       matrix of full measures (nNodes x nMeasures)
%               globalMetrics,      vector of new metrics (1 x nMeasures)
%               hubScore,           vector of top (cumulative measures) hubs
%
%   nodalMetrics (10) = degree, strength, local efficiency, closeness, 
%   betweenness, z-score, participation, eigenvector, pagerank, semi-metricity
%
%   globalMetrics (6) = strength, GC, clustering, global efficiency, 
%   betweenness, semi-metricity
%
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

if nargin == 2
    n = n;
else
    n = 0;
end

nNodes = size(CIJ,1);
nodalMetrics = zeros(nNodes,8); %matrix output
globalMetrics = zeros(1,8); %row vector output
hubScore = zeros(nNodes,1); %column vector of nodes

%Create scaled CIJ for weighted clustering
CIJalt = CIJ; 
CIJalt(1:nNodes+1:end) = 0;
CIJalt = CIJalt/max(CIJalt(:));

L = weight_conversion(CIJ, 'lengths');      %conversion to lengths
d = distance_wei(L);                        %conversion to distance steps
%% Check measures
% myHeavyMeasures - global values

[~, GC] = get_components(CIJ);              %size of giant component
DC = max(GC)/(nNodes-n);                    %disconnected?
C = clustering_coef_wu(CIJalt);             %clustering
Eglob = efficiency_wei(CIJ);                %global efficiency
bc = betweenness_wei(L);                    %betweenness
[~,nSM,SM,~] = computeSemiMetricity(CIJ);   %semi-metricity

% myHeavyMeasures - nodal values
deg = degrees_und(CIJ);                     %degree
S = sum(CIJ,2);                             %strength  
Eloc = efficiency_wei(CIJ,1);               %local efficiency    
cl = 1./(sum(d,2)./(length(d)-1));          %closeness
[Ci, ~] = modularity_louvain_und(CIJ);      %derivation of modularity 
Z = module_degree_zscore(CIJ, Ci, 0);       %modularity z-score
P = participation_coef(CIJ, Ci);            %modularity participation
v = eigenvector_centrality_und(CIJ);        %eigenvector centrality
pr = pagerank_centrality(CIJ, 0.85);        %pagerank centrality
    
% hubCapsHeavy 
[~, I] = sort(deg, 'descend'); %degree centrality
hubScore(I(1:10)) = hubScore(I(1:10)) + 1; %overall
[~, I] = sort(S, 'descend'); %strength
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(cl, 'descend'); %closeness
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(bc, 'descend'); %node betweenness centrality
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(Z, 'descend'); %z-score
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(P, 'descend'); %participation co-efficient
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(v, 'descend'); %eigenvector
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;
[~, I] = sort(pr, 'descend'); %pagerank
hubScore(I(1:10)) = hubScore(I(1:10)) + 1;

%% Parse remaining outputs

globalMetrics = [mean(S) max(GC) DC mean(C) Eglob mean(bc) SM];

nodalMetrics = [deg' S Eloc cl bc Z P v pr nSM]; 

end

