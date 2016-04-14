function [ flags, noFlags, allFlags ] = lesionMaker( CIJ, XYZ, X, radius )
%LESIONMAKER        Makes a lesion in a healthy connectome
%  
%   [flags, noFlags, allFlags] = lesionMaker(CIJ, XYZ, X, radius);
%
%   Takes out parcels based on their radius from a location
%
%   Inputs:     CIJ,        adjacency matrix
%               XYZ,        co-ordinates of all the parcels (n x 3 matrix)
%               X,          tumour centre location
%               radius,     distance from X to flag parcels
%
%   Outputs:    flags,      locations of parcels within radius distance
%               noFlags,    CIJ without flagged parcels
%               allFlags,   CIJ with only flagged parcels
%
% Michael Hart, University of Cambridge, November 2015

%% 1. Calculate 3D Euclidean distance of each nodes XYZ co-ordinate from X

n = size(XYZ,1);
distance = zeros(n,1);

for i = 1:n; %for every parcel
    distance(i,1) = sqrt(((X(1) - XYZ(i,1))^2) + ((X(2) - XYZ(i,2))^2) + ((X(3) - XYZ(i,3))^2));
end

%end up with distances from co-ordinate

%% 2. Flag parcels within radius

flags = distance<radius; %binary vector of parcels less than the radius
not_flagged = distance>radius; %binary vector of parcels outwith radius

%% 3. Remove / isolate flagged parcels from matrix

noFlags = CIJ; allFlags = CIJ;
noFlags(flags,:) = 0; noFlags(:,flags) = 0; %normal matrix without flags
allFlags(not_flagged,:) = 0; allFlags(:,not_flagged) = 0; %matrix with only peri-tumoural flags


