%% Hub Changes Wrapper 
%
% Defines:
% 1. Group specific hubs based on patients with tumours and 
% controls that have been lesioned (0) in a corresponding manner
% 2. Subcategorises into new, missing, same, expected, & predicted hubs
% 3. Extracts measures from individual tumour subjects for each of these 
% overall hub categories & tests differences 
%
% Pipeline
%
% A: establish hubs
%
% 1. Normalise measures
% 2. Compare baselines
% 3. Identify hubs
% 4. Define changes 
%
% B: test hubs
%
% 1. New hubs | overall hubs | individual networks
% 2. Non-hubs | overall hubs | individual networks
%
% C: visualise
%
% 1. Parse to hives
% 2. Parse to BrainNet
%
% Michael Hart, University of Cambridge, March 2016

%% Define inputs & initialise

% measures for each group(non-normalised)
% see individual wrapper for formation
% metrics_tumour; %combined metrics_sub*_tumour.nodalMetrics;
% metrics_lesion; %as above, but for corresponding lesion
% all 116x116 (zeros in lesion +/- control)

nNodes = size(metrics_tumour, 1);
nSubjects = size(metrics_tumour, 3);
nMeasures = size(metrics_tumour, 2);
sMeasures = 5; %2 7 8 9 10 [selected measures for forming Hubs]

%% A: establish hubs

%% 1. Normalise measures

%hubs in tumours 
normTumourMetrics = zeros(nNodes, nMeasures, nSubjects);
for iSubject = 1:nSubjects %individual subjects
    
    for iMeasure = 1:nMeasures %now normalise measures
        normTumourMetrics(:, iMeasure, iSubject) = myNormalNets(metrics_tumour(:, iMeasure, iSubject));
    end %all measures now normalised
        
end %finished all subjects - 116 nodes x 5 subjects matrix

%hubs in lesions
normLesionMetrics = zeros(nNodes, nMeasures, nSubjects);
for iSubject = 1:nSubjects %individual subjects
    
    for iMeasure = 1:nMeasures %now normalise measures
        normLesionMetrics(:, iMeasure, iSubject) = myNormalNets(metrics_lesion(:, iMeasure, iSubject));
    end %all measures now normalised
       
end %finished all subjects - 116 nodes x 5 subjects matrix

%% 2. Compare baselines

%all measures (nb: closeness & SMP aren't valid in zero lesions)
plotmatrix(mean(normTumourMetrics,3))
plotmatrix(mean(normLesionMetrics,3))
plotmatrix(mean(metrics_lesion,3)) %non-normalised

%selected hub measures
plotmatrix(mean(normTumourMetrics(:,[2 7 8 9 10], :),3))
plotmatrix(mean(normLesionMetrics(:,[2 7 8 9 10], :),3))

%% 3. Identify Hubs

HubsT = zeros(nNodes, nSubjects);
HubsMT = zeros(nNodes, sMeasures, nSubjects);
HubsL = zeros(nNodes, nSubjects);
HubsML = zeros(nNodes, sMeasures, nSubjects);
for iSubject = 1:nSubjects
    %in tumours
    normHubs = normTumourMetrics(:, [2 7 8 9 10], iSubject); %select out measures
    Hubs = hubCapsHeavyTwo(normHubs); %calculate hubs
    HubsT(:,iSubject) = Hubs.overall; %look out overall hubs only
    HubsMT(:,:,iSubject) = Hubs.metrics; %hubs per metric
    %in lesions
    normHubs = normLesionMetrics(:, [2 7 8 9 10], iSubject); %select out measures
    Hubs = hubCapsHeavyTwo(normHubs); %calculate hubs
    HubsL(:,iSubject) = Hubs.overall; %look out overall hubs only
    HubsML(:,:,iSubject) = Hubs.metrics; %hubs per metric
end

%% 4. Define hub changes
%hubs in empirical data but not in synthetic lesions
%based on average prevalence of a node being a new hub over subjects
%i.e. 5 means consistently a new hub (max score)
%output is vector nNodes x 1

% new overall hubs
newHubs = HubsT - HubsL;
newHubs = double(newHubs>0); %prevalence; adjust value based on histogram 
%histogram(newHubs)
newHubs = sum(newHubs, 2); %single vector - can use hubViewerOne now
newHubs = double(newHubs>1); %threshold
nnz(newHubs) %24

% hubs not present in tumours but predicted in controls
missingHubs = HubsT - HubsL; %take difference in hub scores
missingHubs = double(missingHubs<0);
%single
missingHubs = sum(missingHubs, 2); %only those not in patients but predicted in controls
missingHubs = double(missingHubs>2); %adjust value based on histogram 
nnz(missingHubs) %23

% same hubs (hub in both HubsT & HubsL)
sameHubs = double(HubsT>0) .* double(HubsL>0);
sameHubs = sum(sameHubs, 2);
sameHubs = sameHubs>1;
nnz(sameHubs) %17

% expected Hubs (those in the empirical data)
expectedHubs = double(HubsT>0);
expectedHubs = sum(expectedHubs,2);
expectedHubs = expectedHubs>2; 
nnz(expectedHubs) %18

% those predicted from simulated lesions
predictedHubs = double(HubsL>0);
predictedHubs = sum(predictedHubs,2);
predictedHubs = predictedHubs>4;
nnz(predictedHubs) %14

%% B: Compare hubs
%
% overall hubs & individual features
%
%% 1. Extract measures
% use normal tumour metrics individually
% use overall group average hubs

newHubsMetrics = normTumourMetrics; %116 * 12 * 5
%newHubsMetrics = metrics_tumour; %not normalised

new_HM = [];
missing_HM = [];
same_HM = [];
expected_HM = [];
predicted_HM = [];

%Extract values
for iSubject = 1:nSubjects
    for iMeasure = 1:nMeasures
        new_HM(:, iMeasure, iSubject) = newHubsMetrics(newHubs==1, iMeasure, iSubject); %do for individual measures
        missing_HM(:, iMeasure, iSubject) = newHubsMetrics(missingHubs==1, iMeasure, iSubject); 
        same_HM(:, iMeasure, iSubject) = newHubsMetrics(sameHubs==1, iMeasure, iSubject);
        expected_HM(:, iMeasure, iSubject) = newHubsMetrics(expectedHubs==1, iMeasure, iSubject);
        predicted_HM(:, iMeasure, iSubject) = newHubsMetrics(predictedHubs==1, iMeasure, iSubject);
    end
end

%% 2. Test differences

[~, p] = ttest2(mean(new_HM,3), mean(missing_HM,3)) %small differences between new & missing
[~, p] = ttest2(mean(new_HM,3), mean(expected_HM,3))
[~, p] = ttest2(mean(new_HM,3), mean(predicted_HM,3)) %key comparison

[~, p] = ttest2(mean(missing_HM,3), mean(expected_HM,3)) %key comparison
[~, p] = ttest2(mean(missing_HM,3), mean(predicted_HM,3))

[~, p] = ttest2(mean(expected_HM,3), mean(predicted_HM,3))

[~, p] = ttest2(mean(new_HM,3), mean(same_HM,3)) %new & same hubs aren't different
[~, p] = ttest2(mean(missing_HM,3), mean(same_HM,3))

mean(mean(edgeDistances(newHubs==1, newHubs==1))) 
mean(mean(edgeDistances(missingHubs==1, missingHubs==1))) 
mean(mean(edgeDistances(sameHubs==1, sameHubs==1))) 
mean(mean(edgeDistances(expectedHubs==1, expectedHubs==1))) 
mean(mean(edgeDistances(predictedHubs==1, predictedHubs==1))) 

[~, p] = ttest2(mean(edgeDistances(expectedHubs==1, expectedHubs==1)), mean(edgeDistances(missingHubs==1,missingHubs==1)))

%% C: visualise
% hives & BrainNet

%% 1. Parse to Hives

Hives(1,:) = mean(mean(new_HM, 3));
Hives(2,:) = mean(mean(predicted_HM, 3));
Hives(3,:) = mean(mean(missing_HM, 3));
Hives(4,:) = mean(mean(expected_HM, 3));

%Hives(:, 12) = 1 - Hives(:, 12); %correction for -ve SMP values (optional)

for i = 1:12
    d3js(:,i) = Hives(:,i) / sum(Hives(:,i)); 
end

%% 2. Visualise in BrainNet

% .node = XYZ Colors Sizes Labels (6)

% new overall hubs
newHubs = HubsT - HubsL;
newHubs = double(newHubs>0); %adjust value based on histogram of prevalences
%single
newHubs = sum(newHubs, 2); %can use hubViewerOne now
newColors = newHubs;
newColors(newColors<2)=0;
newHubs = double(newHubs>1);

missingHubs = HubsT - HubsL; %take difference in hub scores
missingHubs = double(missingHubs<0);
%single
missingHubs = sum(missingHubs, 2); %only those not in patients but predicted in controls
missingColors = missingHubs;
missingColors(missingColors<3) = 0;
missingHubs = double(missingHubs>2); %adjust value based on histogram 

%Define flags
Colors = zeros(116, 1); %empty vector nNodes long
Colors(double(missingHubs>0)==1)=1; %add flags of 1 if next to tumour
Colors(double(newHubs>0)==1)=2; %add flags of 1 if next to tumour

%Define hubs
nodeSizes = zeros(116, 1);
nodeSizes = newColors + missingColors; %iterate over all 12 measures

brainNetFileHubs = [XYZ Colors nodeSizes ones(116,1)]; %then add region names in excel
save('brainNetFileHubs.node', 'brainNetFileHubs', '-ascii');
