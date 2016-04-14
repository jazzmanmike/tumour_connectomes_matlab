function metricChanges = orderedPercolation( CIJ, order )
%ORDEREDPERCOLATION Pre-specified order to failure e.g. comparative group effects
%
%   metricChanges = orderedPercolation(CIJ, order);
%
%   Inputs:     CIJ,                undirected (binary/weighted) connection matrix
%               order,              pre-specified order (focused, random)
%  
%   Outputs:    metricChanges,      matrix of outputs  
%
% Michael Hart, University of Cambridge, January 2016

%% Check & initialise

nNodes = size(CIJ, 1);

%define baselines (for comparisons)
[~, globalBaseMeasures] = percolationMetricsTwo(CIJ);
globalBaseMeasures = repmat(globalBaseMeasures,[nNodes,1]);
gMeasures = size(globalBaseMeasures, 2);

%initialise outcomes
metricChanges = zeros(nNodes, gMeasures); %node x measure

%% Run error loop

nPercolation = 1; %initialised
percolationMatrix = CIJ; %recycled

for iNode = order; %for pre-specified vector
    percolationMatrix(iNode,:) = 0; percolationMatrix(:,iNode) = 0; %take out nodes
    [~, metricChanges(nPercolation,:)] = percolationMetricsTwo(percolationMatrix, nPercolation); %check measures
    nPercolation = nPercolation + 1; %increment counter
end %end per node loop

%% Parse outputs

metricChanges = metricChanges ./ globalBaseMeasures;

end

