%% wrapper: individual hub plasticity

% 1. Create lesions
% 2. Create normalised metrics
% 3. Define hubs 
% 4. Identify hub changes
% 5. Extract metrics
% 6. Parse to group data
% 7. Test for differences

%% Inputs

% Tumour matrices 
% Concordant tumourVectors
% Control matrices

%% 1. Create lesions

%Decide whether to remove [] or zero (0) tumour vector in CIJ +/- T
%Makes a difference for identifying hubs and distance metrics 
%e.g. closeness centrality & semi-metricity

T_lesion = T;
T_lesion(tumourVector,:) = []; T_lesion(:,tumourVector) = [];
CIJ_lesion = CIJ;
CIJ_lesion(tumourVector,:) = []; CIJ_lesion(:,tumourVector) = []; %lesioned control network

%% 2. Create normalised metrics 

%tumour
tumour_metrics = myHeavyMeasures(T_lesion);
nodal_tumour_metrics = tumour_metrics.nodalMetrics;

normal_tumour_metrics = zeros(size(T_lesion,2),12);
for iMeasure = 1:12    
    normal_tumour_metrics(:,iMeasure) = myNormalNets(nodal_tumour_metrics(:,iMeasure));
end %all measures now normalised

%lesion
lesion_metrics = myHeavyMeasures(CIJ_lesion);
nodal_lesion_metrics = lesion_metrics.nodalMetrics;

normal_lesion_metrics = zeros(size(CIJ_lesion,2),12);
for iMeasure = 1:12    
    normal_lesion_metrics(:,iMeasure) = myNormalNets(nodal_lesion_metrics(:,iMeasure));
end %all measures now normalised

%% 3. Define hubs

specific_tumour_metrics = normal_tumour_metrics(:,[2 7 8 9 10]); %select out measures
tumour_hubs = hubCapsHeavyTwo(specific_tumour_metrics); %calculate hubs
HubsT = tumour_hubs.overall;    

specific_lesion_metrics = normal_lesion_metrics(:,[2 7 8 9 10]); %select out measures
lesion_hubs = hubCapsHeavyTwo(specific_lesion_metrics); %calculate hubs
HubsL = lesion_hubs.overall;

%% 4. Identify hub changes

newHubs = double(HubsT>0) - double(HubsL>0); 
newHubs = double(newHubs>0);
missingHubs = double(HubsT>0) - double(HubsL>0);
missingHubs = double(missingHubs<0);
sameHubs = double(HubsT>0)>=1 & double(HubsL>0)>=1;

%% 5. Extract metrics

HubsMetrics = normal_tumour_metrics;

newMetrics005 = [];
missingMetrics005 = [];
sameMetrics005 = [];

for iMeasure = 1:12
    newMetrics005(:,iMeasure) = HubsMetrics(newHubs==1, iMeasure); %do for individual measures
    missingMetrics005(:,iMeasure) = HubsMetrics(missingHubs==1, iMeasure); %do for individual measures
    sameMetrics005(:,iMeasure) = HubsMetrics(sameHubs==1, iMeasure); %do for individual measures
end

%% 6. Parse individual to group data

newMetrics = [newMetrics001; newMetrics002; newMetrics003; ... 
    newMetrics004; newMetrics005];

missingMetrics = [missingMetrics001; missingMetrics002; ... 
    missingMetrics003; missingMetrics004; missingMetrics005];

sameMetrics = [sameMetrics001; sameMetrics002; sameMetrics003; ... 
    sameMetrics004; sameMetrics005];

%or use means (to account for different numbers of hubs per subject)

newMetrics = [mean(newMetrics001); mean(newMetrics002); ... 
    mean(newMetrics003); mean(newMetrics004); mean(newMetrics005)];

missingMetrics = [mean(missingMetrics001); mean(missingMetrics002); ... 
    mean(missingMetrics003); mean(missingMetrics004); ... 
    mean(missingMetrics005)];

sameMetrics = [mean(sameMetrics001); mean(sameMetrics002); ... 
    mean(sameMetrics003); mean(sameMetrics004); mean(sameMetrics005)];

%can 'work out' expected & predicted metrics
expectedMetrics001 = [newMetrics001; sameMetrics001];
expectedMetrics002 = [newMetrics002; sameMetrics002];
expectedMetrics003 = [newMetrics003; sameMetrics003];
expectedMetrics004 = [newMetrics004; sameMetrics004];
expectedMetrics005 = [newMetrics005; sameMetrics005];
expectedMetrics = [expectedMetrics001; expectedMetrics002; ... 
    expectedMetrics003; expectedMetrics004; expectedMetrics005];

predictedMetrics001 = [missingMetrics001; sameMetrics001];
predictedMetrics002 = [missingMetrics002; sameMetrics002];
predictedMetrics003 = [missingMetrics003; sameMetrics003];
predictedMetrics004 = [missingMetrics004; sameMetrics004];
predictedMetrics005 = [missingMetrics005; sameMetrics005];
predictedMetrics = [predictedMetrics001; predictedMetrics002; ... 
    predictedMetrics003; predictedMetrics004; predictedMetrics005];

%or could do same with means (etc.)

%% 7. now test differences

[~, P] = ttest2(newMetrics, missingMetrics)
[~, P] = ttest2(newMetrics, sameMetrics)
[~, P] = ttest2(missingMetrics, sameMetrics)

[~, P] = ttest2(newMetrics, predictedMetrics) %key comparison
[~, P] = ttest2(missingMetrics, expectedMetrics) %key comparison
