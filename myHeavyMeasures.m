function [ Metrics ] = myHeavyMeasures( CIJ )
%MYHEAVYMEASURES Calculates some BCT summary metrics
%
%   MetricsHeavy = myHeavyMeasures(CIJ);
%
%   Total of 26 measures: scalar(9), vector(13), matrix(5)
%
%   Measures stored as a structure to be parsed to e.g. hubCaps
%
%   Inputs:     CIJ,        weighted connectivity matrix/array
%
%   Outputs:    Metrics,    structure with values
%
% Michael Hart, University of Cambridge, December 2015

%% Initialise

nNodes = size(CIJ,1); %number of nodes
nSubjects = size(CIJ,3); %number of subjects in array

%check if binary
if nnz(unique(CIJ))==1
    disp('input is binary - also a binary version available ')
else
    disp('input is weighted')
end

%check if individual matrix or group array
if nSubjects == 1
    disp('input is single group average')
else
    disp('input is an array of individuals')
end

%define outputs
deg = zeros(nNodes,nSubjects) ; EC = zeros(nNodes,nNodes,nSubjects); S = zeros(nNodes,nSubjects);
M0 = zeros(nNodes,nNodes,nSubjects); kden = zeros(1,nSubjects); C = zeros(nNodes,nSubjects); Calt = zeros(nNodes,nSubjects);
T = zeros(1,nSubjects); Eloc = zeros(nNodes,nSubjects); comps = zeros(1,nNodes,nSubjects); Ci = zeros(1,nNodes,nSubjects); 
Q = zeros(1,nSubjects); r = zeros(1,nSubjects); d = zeros(nNodes,nNodes,nSubjects); b = zeros(nNodes,nNodes,nSubjects);
lambda = zeros(1,nSubjects); closeness = zeros(nNodes,nSubjects); efficiency = zeros(1,nSubjects); 
BC = zeros(nNodes,nSubjects); EBC = zeros(nNodes,nNodes,nSubjects); Z = zeros(nNodes,nSubjects); P = zeros(nNodes,nSubjects); 
v = zeros(nNodes,nSubjects); pr = zeros(nNodes,nSubjects);  nodalSM = zeros(nNodes,nSubjects); 
edgeSM = zeros(nNodes,nNodes,nSubjects); globalSM = zeros(1,nSubjects);

%joint_degree = zeros(m+1,m+1,n);

%% Weighted version

% binary = same as binary versions (e.g. converts weights to binary)

for ii = 1:nSubjects
    %similarity
    A = CIJ(:, :, ii); %same if single matrix; individuals if array
    deg(:,ii) = degrees_und(A); %binary
    S(:,ii) = sum(A,2);
    %joint_degree(:,:,ii) = jdegree(A); %binary
    EC(:,:,ii) = edge_nei_overlap_bu(A); %binary
    M0(:,:,ii) = matching_ind_und(A); %binary
    kden(1,ii) = density_und(A); %binary
    %segregation
    C(:,ii) = clustering_coef_wu(A); 
    Ascaled = A/max(A(:)); %scaled weights for clustering
    Ascaled(1:nNodes+1:end) = 0; %zero diagonal
    Calt(:,ii) = clustering_coef_wu(Ascaled);
    T(1,ii) = transitivity_wu(A);
    Eloc(:,ii) = efficiency_wei(A,1);
    comps(1,:,ii) = get_components(A); %binary
    [Ci(1,:,ii), Q(1,ii)] = modularity_louvain_und(A); %binary
    r(1,ii) = assortativity_wei(A,0);
    %integration
    L = weight_conversion(A, 'lengths');
    [d(:,:,ii), b(:,:,ii)] = distance_wei(L); 
    [lambda(1,ii), efficiency(1,ii)] = charpath(d(:,:,ii)); %binary
    %centrality
    closeness(:,ii) = 1./(sum(d(:,:,ii),2)./(length(d(:,:,ii))-1));
    [EBC(:,:,ii), BC(:,ii)] = edge_betweenness_wei(L);
    Z(:,ii) = module_degree_zscore(A,Ci(1,:,ii),0);
    P(:,ii) = participation_coef(A,Ci(1,:,ii));
    v(:,ii) = eigenvector_centrality_und(A);
    pr(:,ii) = pagerank_centrality(A, 0.85);
    [edgeSM(:,:,ii), nodalSM(:,ii), globalSM(1,ii)] = computeSemiMetricity(A);
end %end for loop of individuals


%% Parse outputs

Metrics.degree = deg;
Metrics.strength = S;
%Metrics.jdegree = joint_degree;
Metrics.edge_overlap = EC;
Metrics.matching = M0;
Metrics.density = kden;
Metrics.clustering = C;
Metrics.calt = Calt;
Metrics.transitivity = T;
Metrics.local_efficiency = Eloc;
Metrics.components = comps;
Metrics.modularity = Ci;
Metrics.Qscore = Q;
Metrics.assortativity = r;
Metrics.distance = d;
Metrics.n_steps = b;
Metrics.charpath = lambda;
Metrics.closeness = closeness;
Metrics.global_efficiency = efficiency;
Metrics.node_betweenness = BC;
Metrics.edge_betweenness = EBC;
Metrics.zscore = Z;
Metrics.participation = P;
Metrics.eigenvector = v;
Metrics.pagerank = pr;
Metrics.nodalSM = nodalSM;
Metrics.edgeSM = edgeSM;
Metrics.globalSM = globalSM;

Metrics.globalMetrics = cat(3, kden, mean(Calt), T, Q, r, lambda, efficiency, globalSM); %scalars
Metrics.nodalMetrics = cat(3, deg, S, C, Calt, Eloc, closeness, BC, Z, P, v, ...
    pr, nodalSM); %vectors
Metrics.edgeMetrics = cat(4, b, EBC, M0, edgeSM); %matrices

%Print out the codes to screen
Metrics.globalCode = {'Density'; 'Clustering'; 'Transitivity'; 'Qscore'; ... 
    'Assortativity'; 'lambda'; 'Global_efficiency'; 'Global_SM'}; 

Metrics.nodalCode = {'Degree'; 'Strength'; 'Clustering'; 'altC';  ... 
    'Local_efficiency'; 'Closeness'; 'Betweenness'; 'Zscore'; ...
    'Participation'; 'Eigenvector'; 'Pagerank'; 'semi_metricity'};

Metrics.edgeCode = {'steps'; 'edge_neighbourhood'; 'matching_index'; ...
    'semi_metricity'};

end

