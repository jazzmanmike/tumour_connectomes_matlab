function [ attackScore ] = focusedAttack( CIJ )
%FOCUSEDATTACK Sequential removal of all nodes in order of parameter
%
%   attackScore = focusedAttack(CIJ);
%
%   Inputs:     CIJ,            weighted connectivity matrix
%                 
%   Outputs:    focusedAttack,  matrix of measures (just global)
%
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

nNodes = size(CIJ, 1);

%define baselines
[nodalBaseMeasures, globalBaseMeasures] = percolationMetricsTwo(CIJ); %first output
nMeasures = size(nodalBaseMeasures, 2);
mMeasures = size(globalBaseMeasures, 2);

%initialise outcomes
attackScore = zeros(nNodes, mMeasures, nMeasures); %node x measure x measure

%% Run error loop

for iMeasure = 1:nMeasures %for each measure
    [~, attackOrder] = sort(nodalBaseMeasures(:, iMeasure), 'descend');
    nAttack = 1;
    FA_mat = CIJ; %recycle matrix

    for iNode = attackOrder'; %order of metric
        FA_mat(iNode,:) = 0; FA_mat(:,iNode) = 0; %zero out nodes
        [~, attackChanges] = percolationMetricsTwo(FA_mat, nAttack); %second output
        attackScore(nAttack, :, iMeasure) = attackChanges ./ globalBaseMeasures;
        nAttack = nAttack + 1;
    end %all nodes in measure finished
    
end %all measures checked for all nodes

end

