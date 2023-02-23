if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
end
% add synQuant java path
mydir  = pwd;
addpath(mydir);
Pij = fullfile(mydir, 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile(mydir,'src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
p0 = fullfile(mydir, 'src_synquant/SynQuantVid_v1.2.4.jar');
javaaddpath(p0);%

% data folder
data_folder = '../';
tif_files = dir(fullfile(data_folder, 'crop_embryo_data_500x500x30x40\embryo_TM*.tif'));
% results folder: change to yours
save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
res_folder = fullfile(save_folder,'crop_embryo_data_500x500x30x40');
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
%% synQuant for the 40 frames
z_mat = cell(numel(tif_files), 1);
id_mat = cell(numel(tif_files), 1);
fMaps = cell(numel(tif_files), 1);
q.minIntensity = 0;
for i=1:numel(tif_files)
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));

    sigma = 1;
    sm_im = imgaussfilt3(org_im,sigma);
    % 3D version
    [zMap, synId, fMap] = Synquant4Embryo_Paramater(sm_im, q);
    
    z_mat{i} = single(zMap);
    id_mat{i} = uint16(synId);
    fMaps{i} = fMap;
end
save(fullfile(res_folder, 'synQuant_res.mat'), 'z_mat', 'id_mat','fMaps','-v7.3');
%% refine results from synQuant
eig_res_2d = cell(numel(tif_files), 1);
eig_res_3d = cell(numel(tif_files), 1);
eig_overlay = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    fprintf('cal priCur %d/%d file\n', i, numel(tif_files));
    tic;
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    synId = id_mat{i};
    sigma = 3; % 2d principal curvature smooth scale is smaller
    [eig2d, ~] = principalCv2d(org_im, synId, sigma, fMaps{i});
    sigma = [5,1];
    [eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps{i});

    eig_res_2d{i} = single(eig2d);
    eig_res_3d{i} = single(eig3d);
    eig_overlay{i} = overlay_cl;
    toc;
end
save(fullfile(res_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d',...
    'eig_res_3d','eig_overlay','-v7.3');
%% region refine based on 3d information (infor inside one frame)
refine_res = cell(numel(tif_files), 1);
gcut_res_folder = fullfile(res_folder,'segment_cells');
if ~exist(gcut_res_folder,'dir')
    mkdir(gcut_res_folder);
end
% load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
% load(fullfile(res_folder, 'synQuant_res.mat'));
for i=1:numel(tif_files)
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    synId = id_mat{i};
    eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
    eigAll{1} = eig_res_2d{i};
    eigAll{2} = eig_res_3d{i};
    % pre-frame
    profile on;
    newIdMap = regionWiseAnalysis3dV2(synId, eigAll,...
        org_im, fMaps{i}, gcut_res_folder);
    profile viewer;
    profile off;    
    refine_res{i} = uint32(newIdMap);
end
save(fullfile(res_folder, 'synQuant_refine_res.mat'), 'refine_res','-v7.3');

% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);
