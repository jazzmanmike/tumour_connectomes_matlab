function [ errorScore, errorOrder, errorThreshold ] = randomError( CIJ )
%RANDOMERROR Effects of random and sequential node removal
%
%   [errorScore, errorOrder, errorThreshold] = randomError(CIJ);
%
%   Inputs:     CIJ,             weighted connectivity matrix
%  
%   Outputs:    errorScore,      array with various measure outputs
%               errorOrder,      order of random error percolation   
%               errorThreshold,  mean threshold to disconnect
%               
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

nNodes = size(CIJ, 1);

%define baselines
[~, globalBaseMeasures] = percolationMetricsTwo(CIJ); %second output
globalBaseMeasures = repmat(globalBaseMeasures,[nNodes,1]); %repeat, bigger
nMeasures = size(globalBaseMeasures, 2); %number of outcome measures

%initialise outcomes
errorScore = zeros(nNodes, nMeasures); %node x measure

%% Run error loop

I = randperm(nNodes);
nError = 1;
RE_mat = CIJ; %recycle matrix

for iNode = I; %random permutation order
        RE_mat(iNode,:) = 0; RE_mat(:,iNode) = 0; %set to zeros
        [~, errorScore(nError, :)] = percolationMetricsTwo(RE_mat, nError); %second output
        nError = nError + 1;
end %end per node loop

%% Parse outputs

errorScore = errorScore ./ globalBaseMeasures;
errorOrder = I';
errorThreshold = randomErrorThreshold(CIJ);

end

