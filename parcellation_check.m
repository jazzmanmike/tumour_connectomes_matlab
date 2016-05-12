%Importing parcellation data to matlab & running some quality control checks
%
%Michael Hart, University of Cambridge, 12 May 2016

%% 1. Define FSL

addpath(sprintf('%s/etc/matlab',getenv('FSLDIR'))); %need to be running matlab from command line

%% 2. Import data
%programme is in FreeSurfer directory
%creates a structure base on image header

template_path = ('path/to/parcellation/template');
mask_path = ('path/to/mask/template');

template = load_nifti(template_path);
mask = load_nifti(mask_path);

%% 3. Confirm number of parcels

unique(template.vol) %should equal number of parcels and be integers
unique(mask.vol) %should be binary

%% 4. Lesion template

system(sprintf('fslmaths %s -mul %s lesioned_template', template_path, mask_path)); 
lesioned_template = load_nifti('lesioned_template.nii.gz');

%% 5. Check values

tumourID = unique(lesioned_template.vol);
tumourID(tumourID==0)=[];

%% 6. Check parcel sizes

parcels = load('path/to/ants_n');
hist(parcels)

%rescale to unit interval

parcel_sizes = (parcels - min(parcels)) ./ (max(parcels) - min(parcels)) ;

%% 7. Check XYZ

XYZ = load('/path/to/ants_xyz'); %may need to transpose

XYZ(:,2) = XYZ(:,2)*-1; %may also need to flip directions

%or, if loaded from ANTS

ants_xyz(end) = [];
XYZ = zeros(length(parcel_sizes,3);
for i = 1:3; 
    XYZ(:,i) = ants_xyz(i:3:end); 
end 


%% 7. Check all with BrainNet

% .node = [XYZ Colors Sizes Labels]

Colors = zeros(length(parcels),1);
Colors(tumour_ID) = 1; %flag tumour in red

Sizes = parcel_sizes'; %in unit intervals

BrainNetParcellationFile = [xyz Colors Sizes ones(length(XYZ),1)];
save('BrainNetParcellation.node', 'BrainNetParcellationFile', '-ascii');
BrainNet_MapCfg('BrainMesh', 'BrainNetParcellation.node', 'Cfg.mat', 'Parcellation_check.tiff'); %save image using preset parameters