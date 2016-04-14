%% Tumour Wrapper 
%
% README
%
% Script for basic network analysis of tumour data
%
% Assumes basic pre-processing using semi-metricity scripts e.g. ME-ICA 
% then parcellate (fmri_spt) then wavelet_correlations (SM)
%
% Ultimately concatenates ppp files for a specific wavelet scale
%
% Basic definitions has details of parcellation (AT116)
%
% Calls a number of functions, all organised in /mike_stuff (& sub-folders)
% e.g. fcn_MI_modwt, PersonNetwork, waveletFigures, randomNetworks,
% myHeavyMeasures, measureStats, metricsSides, networkCostMetrics, 
% hubCapsHeavy, hubViewer
%
% Requires some baseline files e.g. P, CIJ, XYZ
%
% Includes;
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
% Michael Hart, University of Cambridge, February 2016

%% A: Quality control

%% A1. Basic definitions

load('~/Documents/MATLAB/Data/Tumours/T001.mat'); %data [matrices, co-ordinates, CIJ]

% Definitions
% ts = time series
% T = tumour network
% CIJ = control network
nNodes = 116; %parcels 
nSubjects = 1; %subjects

regionDescriptions = {'PRE'; 'F20'; 'CHS'; 'CHG'; 'CHG'; 'CHB'; 'CHB'; 
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

% NB: compared with AT116, X & Y vectors have been multiplied by -1

%sub001
%tumourRegions = [PCIN Q O3 P1 P2 AG PQ];
tumourVector = [47 58 63 74 76 80 82];
tumourLocation001 = [24 -64 41]; 

%sub002
tumourRegions = [V1 Q Ling O3 O2 O1 FUSI SMG AG];
tumourVector = [56 58 60 63 65 67 69 78 80];
tumourLocation002 = [37 -63 5];

%sub003
%tumourRegions = V1 Q Ling O3 O2 O1 AG];
%tumourVector = [56 58 60 63 65 67 80];
tumourLocation003 = [29 -69 10];

%sub004
%tumourRegions = [SMA MCIN PCIN POST PQ PCL];
%tumourVector = [30 45 47 71 82 85];
tumourLocation004 = [12 -23 40];

%sub005
%tumourRegions = [R0 PRE POST P2 SMG PUT PAL HES T1 T2];
tumourVector = [27 29 71 76 78 89 91 96 98 102];
tumourLocation005 = [37 -24 25];

CIJ_lesion = CIJ;
CIJ_lesion(tumourVector,:) = []; CIJ_lesion(:,tumourVector) = []; %lesioned control network

%% A2. Wavelet time series analysis

%% Optional: create wavelet correlation matrices

TR=2.42; %TR
lmin=0.02; %high pass filter
lmax=0.20; %highest frequency of your experiment, e.g. if TR=2s you have lmax=1/(2*TR)=1/(2*2)=0.25 Hz (Nyquist frequency)
    
[FCN_MI, WLT] = fcn_MI_MODWT(ts, TR, lmin, lmax);
[PN, p_value] = Person_Network(ts);

T = WLT{2}; %tumour network, wavelet scale 2
T(T<0) = 0; %set negative correlations to 0 - important / TBA

waveletFigures(ts, WLT);

%% A3. Connectivity consistency

% Prevalence matrix
%P = mean(connectivity(:,:,:) > 0, 3); %1 if connections, therefore average number of subjects with a connection

% Outliers
B = double(T>0);
commons = mean(P(B==1)); %low if odd connections
C = B + eye(size(B));
odds = mean(P(C==0)); %high if subject misses connections

%Cost 
kden = density_und(T);

%Degree distribution 
deg = degrees_und(T);

figure; plot(deg, 'o'); %flat line, hmmm

% TODO - make table of basic global matrix stats (commons, odds, density, degrees)
% Could integrate with below

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
MetricsHeavyCIJ = myHeavyMeasures(CIJ);
MetricsHeavyCIJ_lesion = myHeavyMeasures(CIJ_lesion);
MetricsHeavyT = myHeavyMeasures(T);

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

networkCostMetrics(T); %expand to kden cost

%% C. Hubs

% Weighted, control
HubsHeavyCIJ = hubCapsHeavy(MetricsHeavyCIJ); %add specific measures

% Weighted, control
HubsHeavyCIJ_lesion = hubCapsHeavy(MetricsHeavyCIJ_lesion); %add specific measures

% Weighted, tumour
HubsHeavyT = hubCapsHeavy(MetricsHeavyT); %add specific measures

%hubViewer(HubsHeavyT, XYZ);

%parse
hubsCIJ = HubsHeavyCIJ.overall;
hubs_lesion = HubsHeavyCIJ_lesion.overall;
hubsT = HubsHeavyT.overall;

%view distributions
figure; histogram(hubsCIJ);
figure; histogram(hubs_lesion);
figure; histogram(hubsT);

%print out hubs for each measure
hubMetricsT = HubsHeavyT.metrics;
nMeasures = size(hubMetricsT,2); %number of measures

iMetric = zeros(nNodes, nMeasures); 
for ii = 1:nMeasures %actually 8 currently
    [~, iMetric(:,ii)] = sort(heavyHubMetrics(:,ii), 'descend');
end

[regionDescriptions(iMetric(1:10,1)) regionDescriptions(iMetric(1:10,2)) ... 
    regionDescriptions(iMetric(1:10,3)) regionDescriptions(iMetric(1:10,4)) ...
    regionDescriptions(iMetric(1:10,5)) regionDescriptions(iMetric(1:10,6)) ...
    regionDescriptions(iMetric(1:10,7)) regionDescriptions(iMetric(1:10,8))] 

%% D: Comparisons: healthy versus tumour

costT_CIJ = hubsT - hubsCIJ; %subtract values element wise
nHubChangesT_CIJ = nnz(costT_CIJ);
costDifferenceT_CIJ = sqrt(costT_CIJ.^2); %square values then square root

ID = find(costDifferenceT_CIJ>2); %== max(costDifferenceT_CIJ)); %or >2
regionDescriptions{ID}
% compare raw numbers of different nodes, compare normalised sizes of differences

%% E: Visualisation

%% E1. BrainNet

% .node = XYZ Colors Sizes Labels (6)
% .edge = TDL

%Define flags
Colors = zeros(116,1); %empty vector nNodes long
%Colors(tumourVector)=1; %add flags of 1 if next to tumour

%Define hubs
newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = sum(double(newHubs<0), 2); 
newHubs = double(newHubs>2); %group vector, homogenous, n = 14
nodeSizes = newHubs; %iterate over all 12 measures

brainNetFile004 = [XYZ Colors nodeSizes]; %then add region names in excel
save('brainNetFile005.node', 'brainNetFile004', '-ascii');
save('brainNetFile005.edge', 'T', '-ascii', '-tabs'); %optional for nodal sizes only

%% E2. Circos

% .mat.txt = standard formatting of labels and co-ordinates (8 TDL columns)
% i.e. lobe ABBREVIATION color (3) values (3)

% .links.txt = [side(I) name{I} side(J) name{J} weight color];

% NB: vermis removed
% NNB: Side imported from AT116

circosMat = T;
circosMat(vermisID, :) = 0; circosMat(:, vermisID) = 0; %zero out midline/vermis
avgdeg = ((nNodes*(nNodes-1)/2)*0.1)/nNodes; %0.1 = cost
[~, circosMat] = backbone_wu(circosMat, avgdeg); %avgdeg at 10%
%circosMat = threshold_proportional(circosMat, 0.2); %set value for number of edges (N(N-1)/2)
[I, J] = find(double(triu(circosMat)));

fid = fopen('circos.005.links.txt', 'w'); %start a new blank writable file
for i = 1:length(I) %for each row of I (i.e. all connection pairs)
    I_ID = I(i); 
    J_ID = J(i);
    fprintf(fid,'%s\t %s\t %s\t %s\t %f\t %f\n', Side{I_ID}, ...
        regionDescriptions{I_ID}, Side{J_ID}, regionDescriptions{J_ID}, ...
        circosMat(I_ID, J_ID), (Colors(I_ID) + Colors(J_ID)));
end

%% E3. Metric & semi-metric backbones

drawSMBones(T, XYZ); 

%% F: Advanced analysis

%% F1. Lesion toolbox

% Attacks, errors, modular breakdowns, cascades, degeneracy, redundancy

%% F1a. Focused attack
% individual metrics
attackScore = zeros(1, 8, 10, nSubjects);

for iSubject = 1:nSubjects
    attackScore(:,:,:,iSubject) = focusedAttack(CIJ);
end

% group metrics
[~, globalGroupMetrics] = percolationMetrics(CIJ, ones(1,8), 0);
nAttacks = size(globalGroupMetrics, 2);
nMeasures = size(globalGroupMetrics, 2);
groupAttack = zeros(1, nMeasures, nSubjects);

for iSubject = 1:nSubjects
    groupAttack(1,:,iSubject) = orderedPercolation(CIJ, order);
end

AUC = trapz(attackScore);

%% F1b. Hub attack

attackScoreHubs = zeros(nNodes, gMeasures, nSubjects);
attackOrder = zeros(nNodes,nSubjects);

for iSubject = 1:nSubjects
    [attackScoreHubs(:,:,iSubject), attackOrder(:, iSubject)] = hubAttack(connectivity_masked);
end

%% F1c. Random error

% global random removal
order = randperm(1:116); %use for all

for iSubject = 1:nSubjects
    groupError = orderedPercolation(CIJ, order); 
end

% individual random removal
errorScore = zeros(nNodes, 8, nSubjects);
errorOrder = zeros(nNodes,nSubjects);
errorThreshold = zeros(1,nSubjects);

for iSubject = 1:nSubjects
    [errorScore(:,:,iSubject), errorOrder(:,iSubject), errorThreshold(1,iSubject)] = randomError(CIJ);
end

%% F1d. Node removal

%% Module based percolation (group average)

%% Complexity measures (redundancy / degeneracy)


