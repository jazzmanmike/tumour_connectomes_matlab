function hubViewerOne( hubs, XYZ )
%HUBVIEWERONE Plots axial view of nodes proportional to hub status
%   
%   hubViewerOne(hubs, XYZ);
%
%   Inputs: hubs,       single metric vector
%           XYZ,        Euclidean co-ordinates e.g. 
%
% Michael Hart, University of Cambridge, February 2016

%% Define & initialise

nHubs = size(hubs, 2); %number of hubs form above
nNodes = size(hubs, 1);

%% Draw nodes 

figure1 = figure('Name','hub metrics', 'Position', [0 0 1280 700]); %whole page

subplot_1 = subplot(1,2,1,'Parent', figure1);
hold(subplot_1,'on');

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));
nodeSizes = hubs+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


subplot_2 = subplot(1,2,2,'Parent', figure1, 'Position', [0.5 0.15 0.5 0.7]);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

