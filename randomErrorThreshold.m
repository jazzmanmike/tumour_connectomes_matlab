function errorThreshold = randomErrorThreshold( CIJ )
%RANDOMERRORTHRESHOLD Estimates threshold of random error percolation
%   See Network Science, Albert-Laszlo Barabasi
%   
%   errorThreshold = randomErrorThreshold(CIJ);
%
%   Inputs: CIJ,    weighted connectivity matrix
%
%   Outputs: errorThreshold,    percentage of nodes to remove to percolate
%
% Michael Hart, University of Cambridge, January 2016

%% Start short function
degree_dist = sum(CIJ);
k_mean = mean(sum(CIJ));

k2 = (std(degree_dist).^2) / k_mean; %not squared degree

errorThreshold = 1 - (1/((k2/k_mean)-1));

errorThreshold = 1 - errorThreshold; 

end

