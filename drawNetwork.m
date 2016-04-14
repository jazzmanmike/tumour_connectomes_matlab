function drawNetwork( CIJ, XYZ)
%DRAWNETWORK Draws metric & semi-metric backbones
%   Can compare gross topology
%
%   drawNetwork(CIJ, XYZ);
%
%   Inputs: CIJ,        weighted connectivity matrix
%           XYZ,        Euclidean co-ordinates
%
% Michael Hart, University of Cambridge, February 2016

%% Define & initialise

nNodes = size(CIJ, 1);

%% Create backbones
% uses computeSemiMetricity

%1. Compute semi-metricity for each subject
[edgesSM, nodesSM] = computeSemiMetricity(CIJ);

%2. Create prevalence matrix (same as above if single subject)
P_M = edgesSM==1; %metric edges
P_SM = edgesSM>1 & edgesSM<inf; %semi-metric edges

%3. Keep as weights per backbone
edgesMetric = P_M .* CIJ; %gets metric edges if present in > proportion and scale
edgesSemiMetric = P_SM .* CIJ; %both scale edges according to prevalence

%% Plot manually: metric

figure1 = figure('Name','metric & semi-metric edges', 'Position', [1, 1, 1280, 750]);

%%subplot 1
avgdeg = ((nNodes*(nNodes-1)/2)*0.1)/nNodes; %0.1 = cost
[~, edgesMetric] = backbone_wu(edgesMetric, avgdeg); %avgdeg at 10%
figureEdges = nnz(edgesMetric);

Edges=[];
W=[];
avgT = edgesMetric; %average weights of group for line thickness
threshold = min(avgT(avgT~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if edgesMetric(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avgT(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W-threshold; %some sort of normalisation? Minimal value is now 0
W = W*round(1/(1-threshold))*64; %helping for colormap? 64 shades of gray

x1 = XYZ(Edges(:,1),2);
x2 = XYZ(Edges(:,2),2);
y1 = XYZ(Edges(:,1),3);
y2 = XYZ(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

subplot1 = subplot(1,2,1,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot1,'on');

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(0.1+W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = ceil(4 * tiedrank(nodesSM) / length(nodesSM));
hold on
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end

title({'metric edges at 10% cost'});
xlabel(sprintf('%d edges', figureEdges));

%%subplot 2
[~, edgesSemiMetric] = backbone_wu(edgesSemiMetric, avgdeg); %avgdeg at 10%
figureEdges = nnz(edgesSemiMetric);

Edges=[];
W=[];
avgCIJ = edgesSemiMetric; %average weights of group for line thickness
threshold = min(avgCIJ(avgCIJ~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if edgesSemiMetric(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avgCIJ(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W-threshold; %some sort of normalisation? Minimal value is now 0
W = W*round(1/(1-threshold))*64; %helping for colormap? 64 shades of gray

x1 = XYZ(Edges(:,1),2);
x2 = XYZ(Edges(:,2),2);
y1 = XYZ(Edges(:,1),3);
y2 = XYZ(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

subplot2 = subplot(1,2,2,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot2,'on');

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(0.1+W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = ceil(4 * tiedrank(nodesSM) / length(nodesSM));
hold on
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end

title({'semi-metric edges at 10% cost'});
xlabel(sprintf('%d edges', figureEdges));

end

