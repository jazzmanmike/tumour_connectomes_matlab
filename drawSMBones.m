function drawSMBones( CIJ, XYZ)
%DRAWSMBONES Draws metric & semi-metric backbones
%   Can compare gross topology
%
%   drawSMBones(CIJ, XYZ);
%
%   Inputs: CIJ,        weighted connectivity matrix
%           XYZ,        Euclidean co-ordinates
%
% Michael Hart, University of Cambridge, February 2016

%% Define & initialise

nNodes = size(CIJ, 1);
nSubjects = size(CIJ, 3);

%% Create backbones
% uses computeSemiMetricity

%1. Compute semi-metricity for each subject
edgesSM = zeros(nNodes, nNodes, nSubjects);
nodesSM = zeros(nNodes, nSubjects); %for semi-metric percentages
for iSubject = 1:nSubjects
    [edgesSM(:,:,iSubject), nodesSM(:,iSubject)] = computeSemiMetricity(CIJ(:,:,iSubject));
end

%2. Create prevalence matrix (same as above if single subject)
P_M = mean(edgesSM(:,:,:)==1, 3); %metric edges
P_SM = mean(edgesSM(:,:,:)>1 & edgesSM(:,:,:)<inf, 3); %semi-metric edges
nodesSM = mean(nodesSM(:,:), 2, 'omitnan'); %nodal semi-metric percentage

%3. Threshold prevalence matrix & keep as weights
edgesMetric = (P_M > 0.4) .* P_M; %gets metric edges if present in > proportion and scale
edgesSemiMetric = (P_SM > 0.5) .* P_SM; %both scale edges according to prevalence

%% Plot manually: metric

figure1 = figure('Name','metric & semi-metric edges', 'Position', [1, 1, 1280, 750]);

%%subplot 1
subplot1 = subplot(2,4,1,'Parent', figure1);
hold(subplot1,'on');
histogram(P_M, 'BinLimits', [0,1]);
title({'prevalence of metric edges'});

%%subplot 2
edgesMetric = (P_M > 0.6) .* P_M; %gets metric edges if present in > proportion and scale
figureEdges = nnz(edgesMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
threshold = min(avgCIJ(avgCIJ~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if edgesMetric(iEdge, jEdge) ~= 0 %if an edge present
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

subplot2 = subplot(2,4,2,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
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

title({'metric edges at 60% prevalence'});
xlabel(sprintf('%d edges', figureEdges));


%%subplot 3
edgesMetric = (P_M > 0.7) .* P_M; %gets metric edges if present in > proportion and scale
figureEdges = nnz(edgesMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
threshold = min(avgCIJ(avgCIJ~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if edgesMetric(iEdge, jEdge) ~= 0 %if an edge present
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

subplot3 = subplot(2,4,3,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot3,'on');

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

title({'metric edges at 70% prevalence'});
xlabel(sprintf('%d edges', figureEdges));

%%subplot 4
edgesMetric = (P_M > 0.8) .* P_M; %gets metric edges if present in > proportion and scale
figureEdges = nnz(edgesMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
threshold = min(avgCIJ(avgCIJ~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if edgesMetric(iEdge, jEdge) ~= 0 %if an edge present
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

subplot4 = subplot(2,4,4,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot4,'on');

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

title({'metric edges at 80% prevalence'});
xlabel(sprintf('%d edges', figureEdges));

%% Plot manually: semi-metric

%%subplot 5
subplot5 = subplot(2,4,5,'Parent',figure1);
hold(subplot5,'on');
histogram(P_SM, 'BinLimits', [0,1]);
title({'prevalence of semi-metric edges'});

%%subplot 6
edgesSemiMetric = (P_SM > 0.9) .* P_SM; %both scale edges according to prevalence
figureEdges = nnz(edgesSemiMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
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

subplot6 = subplot(2,4,6,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot6,'on');

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

title({'semi-metric edges at 90% prevalence'});
xlabel(sprintf('%d edges', figureEdges));

%%subplot 7
edgesSemiMetric = (P_SM > 0.95) .* P_SM; %both scale edges according to prevalence
figureEdges = nnz(edgesSemiMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
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

subplot7 = subplot(2,4,7,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot7,'on');

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

title({'semi-metric edges at 95% prevalence'});
xlabel(sprintf('%d edges', figureEdges));

%%subplot 8
edgesSemiMetric = (P_SM > 0.97) .* P_SM; %both scale edges according to prevalence
figureEdges = nnz(edgesSemiMetric);

Edges=[];
W=[];
avgCIJ = sum(CIJ(:,:,:), 3) ./ nSubjects; %average weights of group for line thickness
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

subplot8 = subplot(2,4,8,'Parent',figure1, 'XTickLabel', [], 'YTickLabel', []);
hold(subplot8,'on');

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

title({'semi-metric edges at 97% prevalence'});
xlabel(sprintf('%d edges', figureEdges));

%% Plot: need to make square
% uses adjcency_plot_und from BCT

%figure1 = figure('Name','metric & semi-metric edges');

%onlyXY = XYZ(:, [2 3]);

%2D metric
%[x, y] = adjacency_plot_und(edgesMetric, onlyXY);  

%subplot1 = subplot(1,2,1,'Parent',figure1);
%hold(subplot1,'on');
%plot(x,y, 'Parent', subplot1);
%title({'metric 2D'});

%2D semi-metric
%[x, y] = adjacency_plot_und(edgesSemiMetric, onlyXY);  

%subplot2 = subplot(1,2,2,'Parent',figure1);
%hold(subplot2,'on');
%plot(x,y, 'Parent', subplot2);
%title({'semi-metric 2D'});


end

