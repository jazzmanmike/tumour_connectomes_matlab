function parcelViewer( XYZ )
%PARCELVIEWER Plots axial, sagital, coronal views of parcel co-ordinates
%   
%   parcelViewer(XYZ);
%
%   Inputs: XYZ,        Euclidean co-ordinates 
%
% Michael Hart, University of Cambridge, May 2016

%% Define & initialise

nNodes = length(XYZ);
nodes = ones(nNodes,1);

%% Draw nodes 

figure1 = figure('Name','parcellation co-ordinates', 'Position', [1 100 1250 447]); %whole page

subplot_1 = subplot(1,3,1,'Parent', figure1);
hold(subplot_1,'on');

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

title({'axial'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,3),'or','MarkerSize', nodes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

title({'coronal'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_2 = subplot(1,3,3,'Parent', figure1);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

title({'sagital'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end
