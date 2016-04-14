%% Hub Changes Wrapper [original | test bed]
%
% Inputs: patients with tumours and controls
%
% Outputs: differences in hubs between synthetic versus empirical lesions
%
% Take synthetic lesions and compares with empirical data identifying new 
% hubs (not predicted) and missing hubs (predicted but not present).
%
% Algorithm
%
% A: establish hubs
% 1. Normalise measures
% 2. Compare baseline differences
% 3. Mean measure hubs
% 4. Individual measure hubs
%
% B: test hubs
% 1. New hubs | overall hubs | CIJ
% 2. Non-hubs | overall hubs | CIJ
% 3. New hubs | overall hubs | individual networks
% 4. Non-hubs | overall hubs | individual networks
% 5. New hubs | individual hubs | individual networks
% 6. Non-hubs | individual hubs | individual networks
%
% C: evaluate cost
%
% D: visualise | parse to BrainNet
%
% Michael Hart, University of Cambridge, March 2016

%% A: establish hubs

% measures for each group(non-normalised)
% see individual wrapper for formation

%metrics_tumour; %combined metrics_sub*_tumour.nodalMetrics;
%metrics_lesion; %as above, but for corresponding lesion

nNodes = size(metrics_tumour, 1);
nSubjects = size(metrics_tumour, 3);

%% 1. Make new mean measures at group level

%plot tumour
nodalMetricsTumour = squeeze(mean(metrics_tumour,3)); %averaged
plotmatrix(nodalMetricsTumour);

%plot lesion
nodalMetricsLesion = squeeze(mean(metrics_lesion,3)); %averaged
plotmatrix(nodalMetricsLesion); %NB: paths from lesion will be inf

%new selected measures tumour (normalised)
for iMeasure = 1:size(nodalMetricsTumour,2)
    normMeasuresTumour(:,iMeasure) = myNormalNets(nodalMetricsTumour(:,iMeasure));
end

normTumourHubMeasures = normMeasuresTumour(:,[2 7 8 9 10]);

%new selected measures lesion (normalised)
for iMeasure = 1:size(nodalMetricsLesion,2)
    normMeasuresLesion(:,iMeasure) = myNormalNets(nodalMetricsLesion(:,iMeasure));
end

normLesionHubMeasures = normMeasuresLesion(:,[2 7 8 9 10]);

%plot
plotmatrix(normTumourHubMeasures); %better distributions
plotmatrix(normLesionHubMeasures); %note still poor except strength & z-score

%% 2. Compare for global differences at baseline

for iMeasure = 1:size(normTumourHubMeasures,2)
    [~, p] = ttest2(normTumourHubMeasures(:,iMeasure), normLesionHubMeasures(:,iMeasure))
end

[~, p] = ttest2(normTumourHubMeasures, normLesionHubMeasures)
[~, p] = ttest2(nodalMetricsTumour, nodalMetricsLesion)

for iMeasure = 1:size(normTumourHubMeasures,2)
    corr(normTumourHubMeasures(:,iMeasure), normLesionHubMeasures(:,iMeasure))
end

%% 3. Make prevalence matrices of hubs using new normalised mean metrics

%new hub wrapper (only uncorrelated normalised measures)
HubsHeavyTTwo = hubCapsHeavyTwo(normTumourHubMeasures);
HubsHeavyLTwo = hubCapsHeavyTwo(normLesionHubMeasures);

%plot comparisons
hubViewerOne(HubsHeavyTTwo.overall, XYZ) %new, normalised, cummulative
hubViewerFive(HubsHeavyTTwo, XYZ) %selected measures, normalised

hubViewerOne(HubsHeavyLTwo.overall, XYZ) %new, normalised, cummulative
hubViewerFive(HubsHeavyLTwo, XYZ) %selected measures, normalised

%% 4. Define new individual hubs (not using means, normalised selected measures)

%hubs in tumours 
normMetrics = zeros(nNodes, nSubjects);
HubsT = zeros(nNodes, nSubjects);
HubsMT = zeros(nNodes, nMeasures, nSubjects);

for iSubject = 1:nSubjects %individual subjects
    metrics = metrics_tumour(:,:,iSubject); %individual subjects metrics
    
    for iMeasure = 1:size(metrics, 2) %now normalise measures
        normMetrics(:,iMeasure) = myNormalNets(metrics(:,iMeasure));
    end %all measures now normalised
    
    normHubs = normMetrics(:,[2 7 8 9 10]); %select out measures
    Hubs = hubCapsHeavyTwo(normHubs); %calculate hubs
    HubsT(:,iSubject) = Hubs.overall; %look out overall hubs only
    HubsMT(:,:,iSubject) = Hubs.metrics; %hubs per metric
    
end %finished all subjects - 116 nodes x 5 subjects matrix

%hubs in lesions
normMetrics = zeros(nNodes, nSubjects);
HubsL = zeros(nNodes, nSubjects);
HubsML = zeros(nNodes, nMeasures, nSubjects);

for iSubject = 1:nSubjects %individual subjects
    metrics = metrics_lesion(:,:,iSubject); %individual subjects metrics
    
    for iMeasure = 1:size(metrics, 2) %now normalise measures
        normMetrics(:,iMeasure) = myNormalNets(metrics(:,iMeasure));
    end %all measures now normalised
    
    normHubs = normMetrics(:,[2 7 8 9 10]); %select out measures
    Hubs = hubCapsHeavyTwo(normHubs); %calculate hubs
    HubsL(:,iSubject) = Hubs.overall; %look out overall hubs only
    HubsML(:,:,iSubject) = Hubs.metrics; %hubs per metric
    
end %finished all subjects - 116 nodes x 5 subjects matrix

%% 5. Define new hubs 
%hubs in empirical data but not in synthetic lesions
%based on average prevalence of a node being a new hub over subjects
%i.e. 5 means consistently a new hub (max score)

newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = sum(double(newHubs>0), 2); %only those in patients but not predicted in controls
histogram(newHubs); 

newHubs = double(newHubs>2); %adjust value based on histogram of prevalences
nnz(newHubs) %confirm number of hubs
hubViewerOne(newHubs, XYZ); %view the hub locations

%% 6. Define non-hubs
% hubs not present in tumours but predicted in controls

nonHubs = HubsT - HubsL; %take difference in hub scores
nonHubs = sum(double(nonHubs<0), 2); %only those not in patients but predicted in controls
histogram(nonHubs); 

nonHubs = double(nonHubs>3); %adjust value based on histogram 
nnz(nonHubs) %confirm number of hubs
hubViewerOne(nonHubs, XYZ); %view the hub locations

%% B: test hubs
%
% overall hubs & overall features
% overall hubs & individual features
% individual hubs & individual features
%
%% 1. New hubs overall CIJ features

costDifferenceT_CIJ = sqrt((HubsT - HubsL).^2); %square values then square root

metrics = squeeze(MetricsHeavyCIJ.nodalMetrics);
for iMeasure = 1:size(metrics, 2) %now normalise measures
    normMetricsCIJ(:,iMeasure) = myNormalNets(metrics(:, iMeasure));
end %all measures now normalised
    
newHubsMetrics = squeeze(normMetricsCIJ(newHubs==1,:,:)); %do for individual measures
newNonHubsMetrics = squeeze(normMetricsCIJ(newHubs==0,:,:)); %do for individual measures

for i = 1:12; [~, P] = ttest2(newHubsMetrics(:,i), newNonHubsMetrics(:,i)), end %no difference 

edgeDistances = euclideanDistances(XYZ);

mean(mean(edgeDistances(newHubs==1, newHubs==1)))
mean(mean(edgeDistances(newHubs==0, newHubs==0)))

[~, P] = ttest2(edgeDistances(newHubs==1, newHubs==1), edgeDistances(newHubs==0, newHubs==0));

%% 2. Non-hubs overall CIJ features

newNonHubsMetrics = squeeze(normMetricsCIJ(nonHubs==1,:,:)); %do for individual measures
newNonNonHubsMetrics = squeeze(normMetricsCIJ(nonHubs==0,:,:)); %do for individual measures

for i = 1:12; [~, P] = ttest2(newNonHubsMetrics(:,i), newNonNonHubsMetrics(:,i)), end %no difference :(

edgeDistances = euclideanDistances(XYZ);

mean(mean(edgeDistances(nonHubs==1, nonHubs==1)))
mean(mean(edgeDistances(nonHubs==0, nonHubs==0)))

[~, P] = ttest2(edgeDistances(newHubs==1, newHubs==1), edgeDistances(newHubs==0, newHubs==0));

%% 3. New hubs | hubs overall | individual metrics

%i. Define hubs
newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = sum(double(newHubs>0), 2); 
newHubs = double(newHubs>1); %group vector, homogenous, n = 14

%ii. Get invididual measures
nodalMetricsTumour = squeeze(metrics_tumour); %individual measures %unnecessary

%iii. Normalise measures (NB: same results if not normalised)
for iSubject = 1:nSubjects
    for iMeasure = 1:size(nodalMetricsTumour, 2) %now normalise measures
        normMetrics(:,iMeasure,iSubject) = myNormalNets(nodalMetricsTumour(:,iMeasure,iSubject));
    end    
end %individual normalised measures

newHubsMetrics = normMetrics;
newNonHubsMetrics = normMetrics; 

%newHubsMetrics = nodalMetricsTumour;
%newNonHubsMetrics = nodalMetricsTumour; 

%iv. Get individual hub measures
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        newHubsMetrics(newHubs==0, iMeasure, iSubject) = 0; %do for individual measures
        newNonHubsMetrics(newHubs==1, iMeasure, iSubject) = 0; %do for individual measures
    end
end

%v. Test differences
for i = 1:12; 
    [~, P] = ttest2(mean(newHubsMetrics(:,i,:),3), mean(newNonHubsMetrics(:,i,:),3))
end 

mean(mean(edgeDistances(newHubs==1, newHubs==1))) 
mean(mean(edgeDistances(newHubs==0, newHubs==0)))

[~, p] = ttest2(mean(edgeDistances(newHubs==1, newHubs==1), mean(edgeDistances(newHubs==0,newHubs==0))))

%% 4. Non-hubs | hubs overall | individual metrics

%i. define hubs
nonHubs = HubsT - HubsL; %take difference in hub scores
nonHubs = sum(double(nonHubs<0), 2); 
nonHubs = double(nonHubs>3); %group vector, homogenous, n = 14

%ii. define individual tumour measures
nodalMetricsTumour = squeeze(metrics_tumour); %individual measures

%iii. normalise measures
for iSubject = 1:nSubjects
    for iMeasure = 1:size(nodalMetricsTumour, 2) %now normalise measures
        normMetrics(:,iMeasure,iSubject) = myNormalNets(nodalMetricsTumour(:,iMeasure,iSubject));
    end    
end %individual normalised measures

nonHubsMetrics = normMetrics; %no difference if use non-normalised metrics
nonNonHubsMetrics = normMetrics; 

%iv. pull out subject specific hub measures
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        nonHubsMetrics(nonHubs==0, iMeasure, iSubject) = 0; %do for individual measures
        nonNonHubsMetrics(nonHubs==1, iMeasure, iSubject) = 0; %do for individual measures
    end
end

tmp = []
tmp1 = []
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        tmp(:,iMeasure,iSubject) = nonHubsMetrics(nonHubs==0, iMeasure, iSubject); %do for individual measures
        tmp1(:,iMeasure,iSubject) = nonNonHubsMetrics(nonHubs==1, iMeasure, iSubject); %do for individual measures
    end
end

tmp = []
tmp1 = []
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        tmp(:,iMeasure,iSubject) = newHubsMetrics(newHubs==1, iMeasure, iSubject); %do for individual measures
        tmp1(:,iMeasure,iSubject) = newNonHubsMetrics(tmp_hubs==0, iMeasure, iSubject); %do for individual measures
    end
end

tmp = []
tmp1 = []
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        tmp(:,iMeasure) = newHubsMetrics(nonHubs==1, iMeasure, iSubject); %do for individual measures
        tmp1(:,iMeasure) = newNonHubsMetrics(tmp_hubs==1, iMeasure, iSubject); %do for individual measures
    end
end

%v. test differences of measures and distances
for i = 1:12; 
    [~, P] = ttest2(mean(nonHubsMetrics(:,i,:),3), mean(nonNonHubsMetrics(:,i,:),3))
end 

mean(mean(edgeDistances(nonHubs==1, nonHubs==1))) %much longer
mean(mean(edgeDistances(nonHubs==0, nonHubs==0)))

[~, p] = ttest2(mean(edgeDistances(nonHubs==1, nonHubs==1), mean(edgeDistances(nonHubs==0,nonHubs==0))))

sameHubs = double(HubsT>0)>=1 & double(HubsL>0)>=1;

%% 5. New hubs | hubs individual | individual metrics

%i. Define hubs
newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = double(newHubs>1); %116x5, approx. 10 nodes per subject

newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = sum(double(newHubs>0), 2); 
newHubs = double(newHubs>2); %group vector, homogenous, n = 14

%ii. Get individual hub measures
newHubsMetrics = normMetrics;
newNonHubsMetrics = normMetrics; 

tmp = []
tmp1 = []
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        tmp(:,iMeasure) = newHubsMetrics(newHubs==0, iMeasure, iSubject); %do for individual measures
        tmp1(:,iMeasure) = newNonHubsMetrics(newHubs==1, iMeasure, iSubject); %do for individual measures
    end
end

%v. Test differences
for i = 1:12; 
    [~, P] = ttest2(mean(newHubsMetrics(:,i,:),3), mean(newNonHubsMetrics(:,i,:),3))
end 

for i = 1:12; 
    [~, P] = ttest2(mean(newHubsMetrics(:,i,:),3), mean(newNonHubsMetrics(:,i,:),3))
end 

mean(mean(edgeDistances(newHubs==1, newHubs==1))) 
mean(mean(edgeDistances(newHubs==0, newHubs==0)))

[~, p] = ttest2(mean(edgeDistances(newHubs==1, newHubs==1), mean(edgeDistances(newHubs==0,newHubs==0))))

nodeProximity001 = tumourDistance(tumourLocation001, XYZ);
nodeProximity002 = tumourDistance(tumourLocation002, XYZ);
nodeProximity003 = tumourDistance(tumourLocation003, XYZ);
nodeProximity004 = tumourDistance(tumourLocation004, XYZ);
nodeProximity005 = tumourDistance(tumourLocation005, XYZ);

tumourProximity = [nodeProximity001 nodeProximity002 nodeProximity003 nodeProximity004 nodeProximity005];

[~, p] = ttest2(tumourProximity(newHubs==1), tumourProximity(newHubs==0))

%% 6. Non-hubs | hubs individual | individual metrics

%i. define hubs
nonHubs = HubsT - HubsL; %take difference in hub scores
nonHubs = double(nonHubs<-1); %116x5, approx. 25 nodes per subject

%ii. pull out subject specific hub measures
for iSubject = 1:nSubjects
    for iMeasure = 1:12
        nonHubsMetrics(nonHubs(:, iSubject)==0, iMeasure, iSubject) = 0; %do for individual measures
        nonNonHubsMetrics(nonHubs(:, iSubject)==1, iMeasure, iSubject) = 0; %do for individual measures
    end
end

%iii. test differences of measures and distances
for i = 1:12; 
    [~, P] = ttest2(mean(nonHubsMetrics(:,i,:),3), mean(nonNonHubsMetrics(:,i,:),3))
end 

tmp1 = [];
tmp2 = [];
for i = 1:5
    tmp1(:,i) = mean(mean(edgeDistances(nonHubs(:,i)==1, nonHubs(:,i)==1))); %much longer
    tmp2(:,i) = mean(mean(edgeDistances(nonHubs(:,i)==0, nonHubs(:,i)==0)));
end
[~,p] = ttest2(tmp1, tmp2)

[~, p] = ttest2(mean(edgeDistances(nonHubs==1, nonHubs==1), mean(edgeDistances(nonHubs==0,nonHubs==0))))

[~, p] = ttest2(tumourProximity(nonHubs==1), tumourProximity(nonHubs==0))

%% C: Cost function analyses
%
%% 1. Cost differences of different hub measures (5)

costDifferenceT_CIJ = sqrt((HubsMT - HubsML).^2); %square values then square root

totalDifferences = squeeze(sum(costDifferenceT_CIJ)); %5 subjects (rows) 5 measures (columns)

%% BrainNet

% .node = XYZ Colors Sizes Labels (6)
% .edge = TDL

Colors = zeros(116,1); %empty vector nNodes long
%Colors(tumourVector)=1; %add flags of 1 if next to tumour

newHubs = HubsT - HubsL; %take difference in hub scores
newHubs = sum(double(newHubs>0), 2); 
newHubs = double(newHubs>2); %group vector, homogenous, n = 14
newHubs = newHubs .* 3;

nonHubs = HubsT - HubsL; %take difference in hub scores
nonHubs = sum(double(nonHubs<0), 2); 
nonHubs = double(nonHubs>3); %group vector, homogenous, n = 14
nonHubs = nonHubs .* 2;

sharedHubs = HubsT .* HubsL;
sharedHubs = double(sum(sharedHubs,2)>0);
%nodeSizes = ceil(4 * tiedrank(newHubs) / length(newHubs)); %in quartiles
nodeSizes = newHubs + nonHubs + sharedHubs;

Labels = zeros(116,1);

brainNet_NewHubs = [XYZ nodeSizes nodeSizes Labels]; %then add region names in excel
save('brainNet_NewHubs.node', 'brainNet_NewHubs', '-ascii');
