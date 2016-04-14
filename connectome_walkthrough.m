%% Tumour functional connectome walk-through
%
% README
%
% Script for basic functional network analysis of tumour data
%
% Assumes basic pre-processing and parcellation
%
% Basic definitions has details of parcellation (AT116)
%
% Calls a number of functions:
% waveletFigures, randomNetworks, myHeavyMeasures, measureStats,  
% metricsSides, networkCostMetrics, percolationMetrics, costPlots,
% hubCapsHeavy, hubViewer, ...
%
% Requires some baseline files:
% connectivity (control network(s)) & XYZ (co-ordinates) & Side (AT116)
%
% Includes:
%
% A: Quality control

% A1. Basic definitions
% A2. Wavelet time series analysis
% A3. Connectivity consistency
% A4. Generation comparison graphs

% B: Network characterisation
%
% B1. Graph theory measures & checks 
% B2. Cost function analysis of measures

% C: Hubs

% D: Comparison with simulated tumours

% E: visualisation

% E1: BrainNet
% E2: Circos
% E3: backbones

% F: Advanced analysis

% F1. Percolation
% F2. Modular percolation
% F3. Cascades
% F4. Complexity
%
% Michael Hart, University of Cambridge, March 2016

%% A: Quality control

%% A1. Basic definitions

load('~/set/path/to/data/Test_001.mat'); %data [ts, XYZ, CIJ]

% Definitions
% ts = time series
% XYZ = co-ordinates
% connectivity = control network(s)

nNodes = 116; %parcels 

regionDescriptions = {'PRE'; 'F2O'; 'CHS'; 'CHG'; 'CHG'; 'CHB'; 'CHB'; 
'CHT'; 'CHT'; 'CHF'; 'CHF'; 'CVL'; 'F3OP'; 'CVCL'; 'CVCU'; 'CVD'; 'CVT'; 
'CVP'; 'CVU'; 'CVN'; 'F3OP'; 'F3T'; 'F3T'; 'F3O'; 'F3O'; 'RO'; 'RO'; 'SMA'; 
'PRE'; 'SMA'; 'OC'; 'OC'; 'F1M'; 'F1M'; 'F1MO'; 'F1MO'; 'GR'; 'GR'; 'IN'; 
'F1'; 'IN'; 'ACIN'; 'ACIN'; 'MCIN'; 'MCIN'; 'PCIN'; 'PCIN'; 'HIP'; 'HIP'; 
'PHIP'; 'F1'; 'PHIP'; 'AMYG'; 'AMYG'; 'V1'; 'V1'; 'Q'; 'Q'; 'LING'; 'LING'; 
'O3'; 'F1O'; 'O3'; 'O2'; 'O2'; 'O1'; 'O1'; 'FUSI'; 'FUSI'; 'POST'; 'POST'; 
'P1'; 'F1O'; 'P1'; 'P2'; 'P2'; 'SMG'; 'SMG'; 'AG'; 'AG'; 'PQ'; 'PQ'; 'PCL'; 
'F2'; 'PCL'; 'CAU'; 'CAU'; 'PUT'; 'PUT'; 'PAL'; 'PAL'; 'THA'; 'THA'; 'HES'; 
'F2'; 'HES'; 'T1'; 'T1'; 'T1P'; 'T1P'; 'T2'; 'T2'; 'T2P'; 'T2P'; 'T3'; 
'F2O'; 'T3'; 'CHSS'; 'CHSS'; 'CHIS'; 'CHIS'; 'CHCL'; 'CHCL'; 'CHCU'; 
'CHCU'; 'CHS'};

vermisID = [12 14 15 16 17 18 19 20]; %set to 0 for circos

left = [1 4 6 8 10 12:14 19:20 22 24 26 28 31 33 35 37 39 40 42 44 46 ...
48 50 53 55 57 59 61 62 64 66 68 70 72 75 77 79 81 83 84 86 88 90 92 94 ...
97 99 101 103 105 106 108 110 112 114 116];  

right = setdiff(1:116, left);

lvec = zeros(116,1); lvec(left)=1;
rvec = zeros(116,1); rvec(right)=1;

% NB: compared with AT116, X & Y vectors multiplied by -1 (for BrainNet)

%subject 001
tumourRegions = {'PCIN'; 'Q'; 'O3'; 'P1'; 'P2'; 'AG'; 'PQ'};
tumourVector = [47 58 63 74 76 80 82];
tumourLocation = [24 -64 41]; 

CIJ = mean(connectivity, 3); %generate mean connectivity matrix
CIJ_lesion = CIJ;
CIJ_lesion(tumourVector,:) = 0; CIJ_lesion(:,tumourVector) = 0; %lesioned control network

%% A2. Wavelet time series analysis

%intialise outcome
scales = 5; %only do for 5 scales
for iScale = 1:scales
    WLT{iScale} = zeros(nNodes, nNodes);
end

%run modwt & wavelet correlation
for iNode = 1:nNodes %for each node
    WJtX = modwt(ts(:, iNode)); %do a modwt
    for jNode = 1:nNodes %compare with all other nodes
        WJtY = modwt(ts(:, jNode)); %also after modwt
        WT_cor = modwt_wcor(WJtX, WJtY); %form correlation of modwt time series scales
       
        for iLevel = 1:scales %set first 5 scales
            WLT{iLevel}(iNode, jNode) = WT_cor(iLevel); %make a wavelet scale correlation matrix
        end %end per scale 
        
    end %end per nodes comparisons
end %end for all nodes

%define tumour network
T = WLT{2}; %tumour network, wavelet scale 2
T(T<0) = 0; %set negative correlations to 0 - important, TBA
T(eye(nNodes)>0) = 1; %set diagonal to 1 - also important for SMP

waveletFigures(ts, WLT);

% also create control / non-wavelet pearson correlation matrix
pearson_net = corr(ts);
pearson_net(eye(nNodes)>0) = 0;

%% A3. Connectivity consistency

% Prevalence matrix
P = mean(connectivity(:,:,:) > 0, 3); %1 if connections, therefore average number of subjects with a connection

% Outliers
B = double(T>0);
commons = mean(P(B==1)); %low if odd connections
C = B + eye(size(B));
odds = mean(P(C==0)); %high if subject misses connections

% Cost 
kden = density_und(T); %almost 100% connected

%% A4. Generate comparison graphs

mComparisons = 100; %number of comparison graphs

%Weighted
[randomAll, randomWeights] = randomNetworks(T, mComparisons); %very similar if degree = 116

%now show matrices
imagesc(randomAll(:,:,10));
imagesc(randomAll(:,:,10));

%% B: Network characterisation

%% B1. Basic BCT measures

% Weighted
MetricsHeavyCIJ = myHeavyMeasures(CIJ); %controls
MetricsHeavyCIJ_lesion = myHeavyMeasures(CIJ_lesion); %synthetic lesions
MetricsHeavyT = myHeavyMeasures(T); %empirical data

%% B1a. Group measure checks & relationships

% 1. Basic checks (mean, std, histograms, etc)
% NB: global metrics doesn't make sense to do (scalars only)

% 12 nodal metrics (need to sum first)
nodalMetricsT = MetricsHeavyT.nodalMetrics; 
nodalMetricsT(isnan(nodalMetricsT))=0; %NB: remove NaNs - important, TBA
nodalCodes = MetricsHeavyT.nodalCode;

statsNodalT = []; %store structures in cells
for iMetric = 1:size(nodalMetricsT,3); %per metric
    disp(nodalCodes(iMetric));
    statsNodalT(:,iMetric) = measureStats(nodalMetricsT(:,:,iMetric));
end

statsCodes = {'Mean'; 'Standard Deviation'; 'Median'; 'Range'; ...
    '25th Percentile'; '50th Percentile'; '75th Percentile'; ...
    'Semi Interquartile Deviation'; 'Number of outliers'};

%write table
measureStatsTable = array2table(statsNodalT, 'VariableNames', nodalCodes, 'RowNames', statsCodes);

%% B1b. all together

nodalMetricsTPlot = squeeze(nodalMetricsT);
plotmatrix(nodalMetricsTPlot);

%% B1c. left/right comparisons 

metricsSides(nodalMetricsT,left);

%% B2. Cost function analysis

networkCostMetrics(T); %up to 30% cost (where most changes have occured)

%% C. Hubs
HubsHeavyCIJ = hubCapsHeavy(MetricsHeavyCIJ); %control
HubsHeavyCIJ_lesion = hubCapsHeavy(MetricsHeavyCIJ_lesion); %synthetic lesion
HubsHeavyT = hubCapsHeavy(MetricsHeavyT); %empirical data

%view distributions
figure; histogram(HubsHeavyCIJ.overall);
figure; histogram(HubsHeavyCIJ_lesion.overall);
figure; histogram(HubsHeavyT.overall);

%print out hubs for each measure (in empirical data)
hubMetricsT = HubsHeavyT.metrics;
nMeasures = size(hubMetricsT,2); %number of measures

iMetric = zeros(nNodes, nMeasures); 
for ii = 1:nMeasures %actually 8 currently
    [~, iMetric(:,ii)] = sort(hubMetricsT(:,ii), 'descend');
end

[regionDescriptions(iMetric(1:10,1)) regionDescriptions(iMetric(1:10,2)) ... 
    regionDescriptions(iMetric(1:10,3)) regionDescriptions(iMetric(1:10,4)) ...
    regionDescriptions(iMetric(1:10,5)) regionDescriptions(iMetric(1:10,6)) ...
    regionDescriptions(iMetric(1:10,7)) regionDescriptions(iMetric(1:10,8))] 

hubViewerOne(HubsHeavyT.overall, XYZ); %view overall hubs in empirical data
hubViewer(HubsHeavyT, XYZ); %view per metric

%% D: Comparisons: healthy versus tumour

%per measure differences
hubMetricsL = HubsHeavyCIJ_lesion.metrics;
costT_lesion = hubMetricsT - hubMetricsL; %subtract binarised values element wise
costDifference = sum(sqrt(costT_lesion.^2)); %square values then square root

%overall new
hubDifferences = HubsHeavyT.overall - HubsHeavyCIJ_lesion.overall;
newHubs = hubDifferences;
newHubs(newHubs<0) = 0;
regionDescriptions(newHubs>0)
hubViewerOne(newHubs, XYZ)

%overall missing
missingHubs = hubDifferences;
missingHubs(missingHubs>0) = 0;
regionDescriptions(missingHubs<0)
hubViewerOne(missingHubs(*-1), XYZ)

%measures (e.g. for comparisons)
newHubMetrics = nodalMetricsTPlot(newHubs>0, :);
missingHubMetrics = nodalMetricsTPlot(missingHubs<0, :);

%% E: Extended resections

%1. determine distances / plot histograms
nodeProximity = tumourDistance(tumourLocation, XYZ);
histogram(nodeProximity)

%2. establish baseline
baselineHubs = HubsHeavyT.overall;

%3. remove tumour 
resectedNetwork = T;
resectedNetwork(tumourVector,:) = 0; resectedNetwork(:,tumourVector) = 0; %lesioned control network
resectedMeasures = myHeavyMeasures(resectedNetwork);
resectedHubs = hubCapsHeavy(resectedMeasures);
resectionHubs = resectedHubs.overall;

%4. remove margin
extendedNetwork = T;
extendedVector = nodeProximity;
extendedVector(extendedVector<40) = 0; %takes out next 7 parcels closest up to 4cm
extendedVector = logical(double(extendedVector==0));
extendedNetwork(extendedVector,:) = 0; extendedNetwork(:,extendedVector) = 0; %lesioned control network
extendedMeasures = myHeavyMeasures(extendedNetwork);
extendedHubs = hubCapsHeavy(extendedMeasures);
extensionHubs = extendedHubs.overall;

%5. view results
hubViewerExtended(baselineHubs, resectionHubs, extensionHubs, XYZ);

%growing hubs
g1 = baselineHubs - resectionHubs;
g2 = resectionHubs - extensionHubs;
growingHubs = double(g1>0) + double(g2>0);
regionDescriptions(growingHubs==2)
diminishingHubs = double(g1<0) + double(g2<0);
regionDescriptions(diminishingHubs==2)

%% E: Visualisation

%% E1. BrainNet

% .node = XYZ Colors Sizes Labels (6)
% .edge = TDL

%Define flags
Colors = zeros(116,1); %empty vector nNodes long
Colors(tumourVector)=1; %add flags of 1 if next to tumour

%Define hubs
nodeSizes = ceil(4 * tiedrank(newHubs) / length(newHubs));

brainNetFile = [XYZ Colors nodeSizes ones(116,1)]; %then add region names in excel
save('brainNetFile.node', 'brainNetFile', '-ascii');
save('brainNetFile.edge', 'T', '-ascii', '-tabs'); %optional for nodal sizes only

%% E2. Circos

% .mat.txt = standard formatting of labels and co-ordinates (8 TDL columns)
% i.e. lobe ABBREVIATION color (3) values (3)

% .links.txt = [side(I) name{I} side(J) name{J} weight color];

% NB: vermis removed
% NNB: Side imported from AT116

% also: .map.txt = [lobe name R G B z1 z2 z3];

circosMat = T; %define

%or
%circosMat = extendedNetwork;

%or 
%model resection
resect = ones(116,1);
resectionVector = nodeProximity;
resectionVector(resectionVector>55) = 0; %adjust for diameter of resection
resectionVector = logical(double(resectionVector==0));
resect = resectionVector;
circosMat(resect,:) = 0; circosMat(:,resect) = 0;

%format matrix: scale to unit interval & take out vermis
circosMat(vermisID, :) = 0; circosMat(:, vermisID) = 0; %zero out midline/vermis
net_max = max(max(circosMat));
net_min = min(min(circosMat));
circosMat = (circosMat - net_min) ./ (net_max - net_min); 

%Make minimum spanning tree [optional]
avgdeg = ((nNodes*(nNodes-1)/2)*0.1)/nNodes; %0.1 = cost
[~, circosMat] = backbone_wu(circosMat, avgdeg); %avgdeg at 10%

%Parse to file
[I, J] = find(double(triu(circosMat)));
fid = fopen('./circos/raw/circos.links.txt', 'w'); %start a new blank writable file
for i = 1:length(I) %for each row of I (i.e. all connection pairs)
    I_ID = I(i); 
    J_ID = J(i);
    if I_ID ~= J_ID
        fprintf(fid,'%c %s %c %s %g %g\n', Side{I_ID}, ...
        regionDescriptions{I_ID}, Side{J_ID}, regionDescriptions{J_ID}, ...
        (Colors(I_ID) + Colors(J_ID) + 1), circosMat(I_ID, J_ID));
    end %end if statement
end

%Optional, if you just want to run things in perl from matlab
cd circos %need to be in this directory with ./etc
system('./parsemap -map ./raw/map.txt -links ./raw/circos.links.txt')
system('circos -conf etc/circos.conf')

%% E3. Metric & semi-metric backbones

drawNetwork(T, XYZ); 

%% F: Advanced analysis: lesion toolbox

% Attacks, errors, modular breakdowns, cascades, degeneracy, redundancy

%% F1. Focused attack
% individual metrics
attackScore = focusedAttack(T);

AUC = trapz(attackScore(), 1:116);

%% F2. Hub attack
[attackScoreHubs, attackOrder] = hubAttack(T);

AUC = [];

%% F3. Random error
[errorScore, errorOrder, errorThreshold] = randomError(T);

%% F4. Module based percolation (group average)

%% F5. Complexity measures (redundancy / degeneracy)

%% F6. Cascades
