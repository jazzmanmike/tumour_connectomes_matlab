function [ edgeDistances, edgeDistanceCategories ] = euclideanDistances( XYZ, CIJ )
%EUCLIDEANDISTANCES Calculates the distances of each edge
%  
%   [edgeDistances, edgeDistanceCategories] = euclideanDistances(XYZ, CIJ);
%
%   Inputs:     XYZ,                        euclidean co-ordinates (n x 3 matrix)
%               CIJ,                        connectivity matrix (optional)
%
%   Outputs:    edgeDistances,              matrix of distances between nodes
%               edgeDistanceCategories,     short, medium, long
%
% Michael Hart, University of Cambridge, January 2016

%% Initalise & define

nNodes = size(XYZ, 1);
edgeDistances = zeros(nNodes, nNodes);

if nargin == 1
    disp('no matrix supplied - making inverse identity matrix')
    CIJ = ones(nNodes, nNodes);
else
    if nnz(unique(CIJ)) == 1
        disp('matrix is binary')
    else
        disp('matrix is weighted - making binary')
        CIJ = double(CIJ>0);
    end
end


%% Run distance loop

for ii = 1:nNodes 
    for jj = 1:nNodes % for every edge
        edgeDistances(ii,jj) = sqrt(((XYZ(ii,1) - XYZ(jj,1))^2) + ((XYZ(ii,2) - XYZ(jj,2))^2) + ((XYZ(ii,3) - XYZ(jj,3))^2));
    end
end

edgeDistances = edgeDistances .* double(CIJ>0);

%% Split into short / medium / long

lengthVector = squareform(edgeDistances); %upper triangle vector
lowerThird = prctile(lengthVector,((1/3)*100));
upperThird = prctile(lengthVector,((2/3)*100));

shortIndices = squareform(lengthVector<lowerThird);
mediumIndices = squareform(lengthVector>lowerThird & lengthVector<upperThird);
longIndices = squareform(lengthVector>upperThird);

shortEdges = shortIndices;
mediumEdges = mediumIndices;
longEdges = longIndices;

%% Parse outputs

edgeDistanceCategories = [shortEdges + (mediumEdges.*2) + (longEdges.*3)];

end

