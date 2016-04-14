function nodeProximity = tumourDistance( tumourLocation, XYZ )
%TUMOURDISTANCE Calculates the distance of each parcel from the tumour core
%  
%   nodeProximity = tumourDistance(tumourLocation, XYZ);
%
%   Inputs:     tumourLocation,             vector of tumour core XYZ
%               XYZ,                        euclidean co-ordinates (n x 3 matrix)
%      
%
%   Outputs:    nodePromity,                distance vector from tumour core
%
% Michael Hart, University of Cambridge, March 2016

%% Run distance calculator

nNodes = size(XYZ, 1);
nodeProximity = zeros(nNodes, 1);
for ii = 1:nNodes 
    nodeProximity(ii) = sqrt(((tumourLocation(1) - XYZ(ii,1))^2) + ((tumourLocation(2) - XYZ(ii,2))^2) + ((tumourLocation(3) - XYZ(ii,3))^2));
end

end

