function d3js = makeAHive( input1, input2, input3 )
%MAKEAHIVE Parse data to make Hives with d3js
%   
%   d3js = makeAHive(input1, input2, input3);
%
%   Inputs:     inputs[1,2,3],     vectors of measures
%
%   Outputs:    d3js,              txt file to enter in Hives
%
% Michael Hart, University of Cambridge, March 2016

%% Define input matrix

%each input is a column
hives(:, 1) = input1;
hives(:, 2) = input2;
hives(:, 3) = input3;

hives = hives ./ (repmat(sum(hives),3,1)); %scale in proportions

%% Write to file
fid = fopen('d3js','w');

parsed = [0, hives(3,1)+hives(2,1), 1, 1, hives(3,2)+hives(2,2), 1, 0;
    0, hives(3,1), hives(3,1)+hives(2,1), 1, hives(3,2), hives(3,2)+hives(2,2), 1;
    0, 0, hives(3,1), 1, 0, hives(3,2), 2;
    
    1, hives(3,2)+hives(2,2), 1, 2, hives(3,3)+hives(2,3), 1, 0;
    1, hives(3,2), hives(3,2)+hives(2,2), 2, hives(3,3), hives(3,3)+hives(2,3), 1;
    1, 0, hives(3,2), 2, 0, hives(3,3), 2;
    
    2, hives(3,3)+hives(2,3), 1, 0, hives(3,1)+hives(2,1), 1, 0;
    2, hives(3,3), hives(3,3)+hives(2,3), 0, hives(3,1), hives(3,1)+hives(2,1), 1
    2, 0, hives(3,3), 0, 0, hives(3,1), 2];

fprintf(fid, '{source: {x: %g, y0: %g, y1:  %g}, target: {x: %g, y0: %g, y1:  %g}, group:  %g},\n', parsed');

d3js = fid;
end

