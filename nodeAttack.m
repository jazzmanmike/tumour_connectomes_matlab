function [ deltaNode, deltaGlobal ] = nodeAttack( CIJ )
%NODEATTACK Effects of single node removal
%
%   [deltaNode, deltaGlobal] = nodeAttack(CIJ);
%
%   Inputs:     CIJ,            weighted connectivity matrix
%  
%   Outputs:    deltaNode,      percentage change in nodal measures 
%               deltalGlobal,   change in global measures
%
%   NB: buglet - column 3 of deltaNode will always be 0
%
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

%define constants
nNodes = size(CIJ, 1); %number of parcels/nodes

%define baselines
[nodalBaseMeasures, globalBaseMeasures] = percolationMetricsTwo(CIJ); %second output
%globalBaseMeasures = repmat(globalBaseMeasures,[nNodes,1]);

nMeasures = size(nodalBaseMeasures, 2);
mMeasures = size(globalBaseMeasures, 2);

%initialise outcome
deltaNode = zeros(nNodes, nMeasures, nNodes);
deltaGlobal = zeros(nNodes, mMeasures); 

%% Run attacks

for iNode = 1:nNodes; %per node

    W = CIJ; %start each time with fresh matrix
    W(iNode,:) = 0; W(:,iNode) = 0; %set all edges of node to 0
    [nodeChange, globalChange] = percolationMetricsTwo(W, 1); %calculate new nodal & global measures
    
    deltaNode(:,:,iNode) = ((nodalBaseMeasures - nodeChange) ./ nodalBaseMeasures) .* 100; %convert to percentage changes
    deltaGlobal(iNode,:) = ((globalBaseMeasures - globalChange) ./ globalBaseMeasures) .* 100;

end %end node loop

%% Parse output

%deltaGlobal = ((globalBaseMeasures - deltaGlobal) ./ globalBaseMeasures) .* 100; 
%deltaNode = ((nodalBaseMeasures - deltaNode) ./ nodalBaseMeasures) .* 100;

%NB MSM = 0 if both metric, NaN if both semi-metric, 1 if become semi-metric, -Inf if become metric
end

