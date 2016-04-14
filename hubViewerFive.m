function hubViewerFive( hubs, XYZ )
%HUBVIEWERFIVE Plots axial view of nodes proportional to hub status
%   
%   hubViewerFive(hubs, XYZ);
%
%   Inputs: hubs,       structure from myHeavyHubs (selects .metrics)
%           XYZ,        Euclidean co-ordinates e.g. 
%
% Michael Hart, University of Cambridge, February 2016

%% Define & initialise
hubs = hubs.metrics; %parse from string
nHubs = size(hubs, 2); %number of hubs form above
nNodes = size(hubs, 1);
hubNames = {'strength'; 'betweenness'; 'zscore'; 'participation'; 'eigenvector'};

%% Draw nodes 

figure1 = figure('Name','hub metrics', 'Position', [100 100 1280 700]); %whole page

for iHub = 1:nHubs %for each of 5 hubs, make a subplot
    subplot_1_{iHub} = subplot(2,5,iHub,'Parent', figure1);
    hold(subplot_1_{iHub},'on');
    nodeSizes = ceil(4 * tiedrank(hubs(:,iHub)) / length(hubs));
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
    end %end individual hub loop

    set(gca,'xaxislocation','top');
    title(hubNames{iHub});
    xlabel({'anterior'});
    ylabel({'lateral'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');
end

for iHub = 1:nHubs
    subplot_{iHub+5} = subplot(2,5,iHub+5,'Parent', figure1);
    hold(subplot_{iHub+5},'on');
    nodeSizes = ceil(4 * tiedrank(hubs(:,iHub)) / length(hubs));
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
    end %end individual hub loop

    set(gca,'xaxislocation','top');
    title(hubNames{iHub});
    xlabel({'superior'});
    ylabel({'posterior'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all hubs

end

