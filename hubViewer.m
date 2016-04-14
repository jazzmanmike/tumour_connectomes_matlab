function hubViewer( hubs, XYZ )
%HUBVIEWER Plots axial view of nodes proportional to hub status
%   
%   hubViewer(hubs, XYZ);
%
%   Inputs: hubs,       structure from myHeavyHubs (selects .metrics)
%           XYZ,        Euclidean co-ordinates e.g. 
%
% Michael Hart, University of Cambridge, February 2016

%% Define & initialise
hubs = hubs.metrics; %parse from string
nHubs = size(hubs, 2); %number of hubs form above
nNodes = size(hubs, 1);
hubNames = {'degree'; 'strength'; 'closeness'; 'betweenness'; 'zscore'; 'participation'; 'eigenvector'; 'pagerank'};

%% Draw nodes 

figure1 = figure('Name','hub metrics', 'Position', [100 100 1280 700]); %whole page

for iHub = 1:nHubs %for each of 8 hubs, make a subplot
    subplot_{iHub} = subplot(2,4,iHub,'Parent', figure1);
    hold(subplot_{iHub},'on');
    nodeSizes = ceil(4 * tiedrank(hubs(:,iHub)) / length(hubs));
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*4, 'MarkerEdgeColor','k','MarkerFaceColor','r');
    end %end individual hub loop

    set(gca,'xaxislocation','top');
    title(hubNames{iHub});
    xlabel({'anterior'});
    ylabel({'lateral'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all hubs

end

