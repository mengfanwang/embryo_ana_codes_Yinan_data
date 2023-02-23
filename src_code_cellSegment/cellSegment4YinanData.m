if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
    save_folder = '/home/congchao/Desktop/embryo_res_folder';
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
    save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
end

addpath('src_code_matlab');
addpath('src_code_cellSegment');
addpath('src_code_cellTracker');


%mydir  = '/home/ccw/Desktop/yinan_data_samples/time_000';
mydir  = '/home/congchao/storage/Yinan_data/TM0-49/crop';
addpath('./');
% data folder
data_folder = mydir;
tif_files = dir(fullfile(data_folder, '/*.tif'));
% results folder
res_folder = fullfile(save_folder,'yinan_data_res');
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
minIntensity = 500; % The middle of two Gaussian intensity distributions (
                    % manually set but should learn from data)

%% synQuant
% add synQuant java path
Pij = fullfile('src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
%p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.4.jar');
p0 = fullfile('src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);%

z_mat = cell(numel(tif_files), 1);
id_mat = cell(numel(tif_files), 1);
fMaps = cell(numel(tif_files), 1);
q.minIntensity = minIntensity;
for i=1:numel(tif_files)
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    tic;
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [h, w, slices] = size(org_im);
    %out_ims = SliceImage(in_im);
    %org_im = imresize3(org_im,round([h/2 w/2 slices/2]));
    sigma = 1;
    sm_im = imgaussfilt3(org_im,sigma);
    %q.posEigMap = eig_res_3d{i}>0;
    % 3D version
    [zMap, synId, fMap] = Synquant4Embryo_Paramater(sm_im, q);
    
    z_mat{i} = single(zMap);
    id_mat{i} = uint16(synId);
    fMaps{i} = fMap;
    toc;
end
save(fullfile(res_folder, 'synQuant_res.mat'), 'z_mat', 'id_mat','fMaps','-v7.3');
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);

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
    fMaps = org_im > minIntensity;
    [eig2d, ~] = principalCv2d(org_im, synId, sigma, fMaps);
    sigma = [5,1];
    [eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps);

    eig_res_2d{i} = single(eig2d);
    eig_res_3d{i} = single(eig3d);
    eig_overlay{i} = overlay_cl;
    toc;
end
save(fullfile(res_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d',...
    'eig_res_3d','eig_overlay','-v7.3');

%% calculate the variance map of all frames
varMap = cell(numel(tif_files), 1);
scale_term = 1100;
for i=1:numel(tif_files)
    disp(i);
    tic;
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    vid = 255*org_im/scale_term;
    varMap{i} = cell(3,2);
    [varMap{i}{1,1}, varMap{i}{2,1},varMap{i}{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8, 3);
    vid_stb = sqrt(vid+3/8);
    [varMap{i}{1,2}, varMap{i}{2,2},varMap{i}{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8, 3);
    toc;
end
save(fullfile(res_folder, 'varianceMap.mat'), 'varMap','-v7.3');

%% region refine based on 4d information (infor across >1 frame)
load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(res_folder, 'synQuant_res.mat'));
load(fullfile(res_folder, 'varianceMap.mat'));
refine_res = cell(numel(tif_files), 1);
threshold_res = cell(numel(tif_files), 1);
multi_frames_flag = false; % use multiple frames for segmentation
cell_wise_save_flag = false; % save each cell segment
for i=1:numel(tif_files)
    tic;
    if cell_wise_save_flag
        seg_res_folder = fullfile(res_folder,'segment_cells');
        seg_res_folder = fullfile(seg_res_folder,['frame_', num2str(i)]);
        if ~exist(seg_res_folder,'dir')
            mkdir(seg_res_folder);
        end
    end
    if multi_frames_flag % use 3 consecutive frames
        org_im = cell(3,1);
        synId = cell(3,1);
        eigAll = cell(3,1);
        varMapAll = cell(3,1);
        for j=i-1:i+1
            if j>0 && j<=numel(tif_files)
                tmp_im = tifread(fullfile(tif_files(j).folder, tif_files(j).name));
                tmp_im = 255*tmp_im/scale_term;
                org_im{j-i+2} = tmp_im;
                synId{j-i+2} = id_mat{j};
                eigAll{j-i+2} = cell(2,1); % save both 2d and 3d pincipal curvature
                eigAll{j-i+2}{1} = eig_res_2d{i};
                eigAll{j-i+2}{2} = eig_res_3d{i};
                varMapAll{j-i+2} = varMap{j};
            end
        end
    else
        tmp_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
        org_im = 255*tmp_im/scale_term;
        synId = id_mat{i};
        eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
        eigAll{1} = eig_res_2d{i};
        eigAll{2} = eig_res_3d{i};
        varMapAll = varMap{i};
    end
    
    %profile on;
    [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
            eigAll,org_im, varMapAll, []);%, i
    %profile viewer;
    %profile off;
    toc;
    refine_res{i} = uint32(newIdMap);
    threshold_res{i} = uint8(thresholdMap);
end
save(fullfile(res_folder, 'synQuant_refine_res_4d_v9.mat'), 'refine_res',...
    'threshold_res','-v7.3');
% v2: grow fg if foreground touching boundary in xy direction
% v3: grow fg if the detected cell touches boundary in xy direction
% v4: grow fg if the detected cell touches boundary in any 3 direction
% v5: grow fg in corresponding direction if the detected cell touches x, y
%     or z direction
% v6: shift = [10 10 2] ==> [20 20 4]
% v7: test only the connected component with one seed, remove other seed
% regions (this indeed does not need 3 frames)
% v8: rank all the seed region based on their intensity, remove brighter
% regions before refine dimmer ones
% v9: use 2d principal curvautre to re-detect seeds in splitFGintoCells.m
%% second time of synQuant
% add synQuant java path
Pij = fullfile('../../ParticleTracking', 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('../../ParticleTracking','src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
%p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.4.jar');
p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);%
load(fullfile(res_folder, 'synQuant_refine_res_4d_v9.mat'), 'refine_res');
id_mat_2nd = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    tic;
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    sigma = 1;
    sm_im = imgaussfilt3(org_im,sigma);
    posEigMap = eig_res_3d{i}>0;
    % 3D version
    id_mat_2nd{i} = Synquant4Embryo_2iter(sm_im, refine_res{i}, posEigMap);
    toc;
end
save(fullfile(res_folder, 'synQuant_res_2iter_indirect4v9.mat'), 'id_mat_2nd','-v7.3');
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);

%% refine the new detected cells
load(fullfile(res_folder, 'synQuant_res_2iter_indirect4v9.mat'), 'id_mat_2nd');
load(fullfile(res_folder, 'synQuant_refine_res_4d_v9.mat'));
seg_res_folder = [];
id_mat_com = refine_res;
for i=1:numel(tif_files)
    cell_append = max(id_mat_2nd{i}(:));
    valid_fg = id_mat_com{i}>0;
    id_mat_com{i}(valid_fg) = id_mat_com{i}(valid_fg) + cell_append;
    id_mat_com{i} = id_mat_com{i} + id_mat_2nd{i};
end
refine_res_1st = refine_res;
threshold_res_1st = threshold_res;
%svf = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\segmentation_from_synQuant\';
multi_frames_flag = false; % use multiple frames for segmentation
for i=1:numel(tif_files)
    if multi_frames_flag % use 3 consecutive frames
        org_im = cell(3,1);
        synId = cell(3,1);
        eigAll = cell(3,1);
        varMapAll = cell(3,1);
        for j=i-1:i+1
            if j>0 && j<=numel(tif_files)
                tmp_im = tifread(fullfile(tif_files(j).folder, tif_files(j).name));
                tmp_im = 255*tmp_im/scale_term;
                org_im{j-i+2} = tmp_im;
                synId{j-i+2} = id_mat_com{j};
                eigAll{j-i+2} = cell(2,1); % save both 2d and 3d pincipal curvature
                eigAll{j-i+2}{1} = eig_res_2d{i};
                eigAll{j-i+2}{2} = eig_res_3d{i};
                varMapAll{j-i+2} = varMap{j};
            end
        end
    else
        tmp_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
        org_im = 255*tmp_im/scale_term;
        synId = id_mat_com{i};
        eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
        eigAll{1} = eig_res_2d{i};
        eigAll{2} = eig_res_3d{i};
        varMapAll = varMap{i};
    end
    
    tic;
    test_ids = 1:double(max(id_mat_2nd{i}(:)));
    [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
            eigAll,org_im, varMapAll, test_ids, []);
    toc;
    
    valid_fg = refine_res_1st{i}>0;
    refine_res_1st{i}(valid_fg) = refine_res_1st{i}(valid_fg) + ...
        uint32(max(newIdMap(:)));
    refine_res{i} = uint32(newIdMap) + refine_res_1st{i};
    threshold_res{i} = uint8(thresholdMap) + threshold_res_1st{i};
end
save(fullfile(res_folder, 'synQuant_refine_res_4d_v9plus.mat'), 'refine_res',...
    'threshold_res','-v7.3');

% v6+1: repeat one more time to find missed cells
