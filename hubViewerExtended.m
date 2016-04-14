function hubViewerExtended( hubs1, hubs2, hubs3, XYZ )
%HUBVIEWERONE Plots axial view of nodes proportional to hub status
%   
%   hubViewerOne(hubs1, hubs2, hubs3, XYZ);
%
%   Inputs: hubs1,      empirical data
%           hubs2,      tumour resected
%           hubs3,      extended resection
%
%           XYZ,        Euclidean co-ordinates e.g. 
%
% Michael Hart, University of Cambridge, March 2016

%% Define & initialise
nNodes = size(hubs1, 1);

%% First hubs 

figure1 = figure('Name','hub metrics', 'Position', [100 100 1240 700]); %whole page

subplot_1 = subplot(2,3,1,'Parent', figure1);
hold(subplot_1,'on');

%nodeSizes = ceil(4 * tiedrank(hubs1) / length(hubs1));
nodeSizes = hubs1+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_2 = subplot(2,3,4,'Parent', figure1);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%% Second hubs 

subplot_3 = subplot(2,3,2,'Parent', figure1);
hold(subplot_3,'on');

%nodeSizes = ceil(4 * tiedrank(hubs2) / length(hubs2));
nodeSizes = hubs2+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_4 = subplot(2,3,5,'Parent', figure1);
hold(subplot_4,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%% Third hubs 

subplot_5 = subplot(2,3,3,'Parent', figure1);
hold(subplot_5,'on');

%nodeSizes = ceil(4 * tiedrank(hubs3) / length(hubs3));
nodeSizes = hubs3+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_6 = subplot(2,3,6,'Parent', figure1);
hold(subplot_6,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

