function [ attackScoreHubs, attackOrder ] = hubAttack( CIJ )
%HUBATTACK Focused attack based on hubs
%
%   [attackScoreHubs, attackOrder] = hubAttack(CIJ);
%
%   Inputs:     CIJ,                weighted connectivity matrix
%                 
%   Outputs:    attackScoreHubs,    matrix of measures (nNodes x nMeasures)
%               attackOrder,        order of percolation (based on hubs)
%              
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

nNodes = size(CIJ, 1);

%define baselines
[~, globalBaseMeasures, hubScore] = percolationMetrics(CIJ); %first output
globalBaseMeasures = repmat(globalBaseMeasures, nNodes, 1);
gMeasures = size(globalBaseMeasures, 2);

%initialise outcomes
attackScoreHubs = zeros(nNodes, gMeasures); %node x measure 

%% Run attacks

I = find(hubScore==max(hubScore)); %index of a top hub
if length(I)>1
    I = datasample(I,1); %randomly selects a top hub if >1
end

attackOrder = I;
nAttack = 1;
FA_mat = CIJ; %recycle matrix

for iNode = 1:nNodes; %iterations for all nodes
    FA_mat(I,:) = 0; FA_mat(:,I) = 0; %make zeros for hub node
    [~, attackScoreHubs(iNode, :), hubScore] = percolationMetrics(FA_mat, nAttack); %second output
    I = find(hubScore==max(hubScore)); 
    
    if length(I)>1 %see above commentary for if loop
        I = datasample(I,1); 
    end
    
    if length(I)==0 %in case no hubs left
        remainingNodes = (1:116)'; %column of scalars 1:116
        remainingNodes(attackOrder)=0; %zero out those nodes already gone
        I = datasample(remainingNodes,1);
    end
    
    nAttack = nAttack + 1;
    attackOrder = [attackOrder; I];
end %end per node loop

%% Parse outputs

attackScoreHubs = attackScoreHubs ./ globalBaseMeasures;

end

