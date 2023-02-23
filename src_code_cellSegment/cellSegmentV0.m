if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
end
addpath('src_code_matlab');
addpath('src_code_cellSegment');
addpath('src_code_cellTracker');

addpath('../../ParticleTracking/src_code/');
addpath('../../ParticleTracking/src_uTrack/');
addpath('../../ParticleTracking/CINDA/');

mydir  = '';
addpath(mydir);
% data folder
data_folder = mydir;
tif_files = dir(fullfile(data_folder, 'crop_embryo_data_500x500x30x40\embryo_TM*.tif'));
% results folder
save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
res_folder = fullfile(save_folder,'crop_embryo_data_500x500x30x40');
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end

%% synQuant for the 40 frames
% add synQuant java path
Pij = fullfile(mydir, 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile(mydir,'src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
p0 = fullfile(mydir, 'src_synquant/SynQuantVid_v1.2.4.jar');
javaaddpath(p0);%

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
%% calculate the variance map of all frames
varMap = cell(numel(tif_files), 1);
scale_term = 5000;
parfor i=1:numel(tif_files)
    disp(i);
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    vid = 255*org_im/scale_term;
    varMap{i} = cell(3,2);
    [varMap{i}{1,1}, varMap{i}{2,1},varMap{i}{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8);
    vid_stb = sqrt(vid+3/8);
    [varMap{i}{1,2}, varMap{i}{2,2},varMap{i}{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8);
end
save(fullfile(res_folder, 'varianceMap.mat'), 'varMap','-v7.3');

%% region refine based on 3d information (infor inside one frame)
refine_res = cell(numel(tif_files), 1);
load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(res_folder, 'synQuant_res.mat'));
scale_term = 5000;
profile on;
for i=1:numel(tif_files)
    seg_res_folder = fullfile(res_folder,'segment_cells');
    seg_res_folder = fullfile(seg_res_folder,['frame_', num2str(i)]);
    if ~exist(seg_res_folder,'dir')
        mkdir(seg_res_folder);
    end
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    org_im = 255*org_im/scale_term;
    %nxt_im = tifread(fullfile(tif_files(i+1).folder, tif_files(i+1).name));
    synId = id_mat{i};
    eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
    eigAll{1} = eig_res_2d{i};
    eigAll{2} = eig_res_3d{i};
    % pre-frame
    newIdMap = regionWiseAnalysis3dV2(synId, eigAll,...
        org_im, fMaps{i}, varMap{i}, seg_res_folder);
%     newIdMap = regionWiseAnalysis3dV2(synId, eigAll,...
%         org_im, fMaps{i}, varMap{i}, []);

    refine_res{i} = uint32(newIdMap);
end
profile viewer;
profile off;
save(fullfile(res_folder, 'synQuant_refine_res.mat'), 'refine_res','-v7.3');

%% The previous sections are for cell segmentation; the following sections are not used %%
% -----------The tracking part, please refer to track_refine.m-----------------------%
% ----------------------------------------------------------May.19th,2020------------%

%% region refine based on 4d information (infor across >1 frame)
load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(res_folder, 'synQuant_res.mat'));
load(fullfile(res_folder, 'varianceMap.mat'));
scale_term = 5000;
refine_res = cell(numel(tif_files), 1);
threshold_res = cell(numel(tif_files), 1);
seg_res_folder = [];
for i=1:numel(tif_files)
    disp(i);
%     seg_res_folder = fullfile(res_folder,'segment_cells');
%     seg_res_folder = fullfile(seg_res_folder,['frame_', num2str(i)]);
%     if ~exist(seg_res_folder,'dir')
%         mkdir(seg_res_folder);
%     end
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
    % pre-frame
    tic;
%     [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
%         eigAll,org_im, varMapAll, seg_res_folder);
    [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
        eigAll,org_im, varMapAll, []);
    toc;
    refine_res{i} = uint32(newIdMap);
    threshold_res{i} = uint8(thresholdMap);
end
save(fullfile(res_folder, 'synQuant_refine_res_4d.mat'), ...
    'refine_res','threshold_res','-v7.3');

[h,w,~] = size(threshold_res{1});
outThresMap = uint8(zeros(h,w,numel(threshold_res)));
vidMaxProj = uint8(zeros(h,w,numel(threshold_res)));
for i=1:numel(threshold_res)
    tmp_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    tmp_im = uint8(255*tmp_im/scale_term);
    outThresMap(:,:,i) = max(threshold_res{i},[], 3);
    vidMaxProj(:,:,i) = max(tmp_im,[], 3);
end
tifwrite(vidMaxProj, fullfile(res_folder, 'vidMapProj'));
tifwrite(outThresMap, fullfile(res_folder, 'threMap'));

%% tracking way 2: considering overlapping ratios
if ~exist('refine_res','var')
    load(fullfile(res_folder, 'synQuant_refine_res.mat'));
end
% let's first downsample the detection results
if ~exist('org_refine_res', 'var')
    org_refine_res = refine_res;
else
    refine_res = org_refine_res;
end
[h, w, z] = size(refine_res{1});
ds_sz = round([h/3, w/3, z]);
for i = 1:numel(refine_res)
    refine_res{i} = imresize3(refine_res{i},ds_sz,'method','nearest');
end
% main function
[movieInfo,movieInfoAll] = mcfTracking_cell(refine_res);
movieInfo.xCoord = movieInfo.orgCoord(:,1);
movieInfo.yCoord = movieInfo.orgCoord(:,2);
movieInfo.zCoord = movieInfo.orgCoord(:,3);
save(fullfile(res_folder, 'movieInfo.mat'),'movieInfo','movieInfoAll');


%% validate the results
if ~exist('movieInfo','var')
    load(fullfile(res_folder, 'movieInfo.mat'),'movieInfo');
end
load('crop_embryo_data_500x500x30x40\gt.mat','gt_mat');
gt_tracks = loc2id(gt_mat, org_refine_res);
[forward_mat, ratio_mat, track_head] = validate_res(movieInfo, gt_tracks);
zz = round(nanmax(ratio_mat,[],2),2);
%% visualize the reuslts
disp_folder = 'synQuant_refine_res';
max_cell_num = max(cellfun(@(x) max(x(:)), refine_res));
cmap = jet(double(max_cell_num));
cmap = cmap(randperm(max_cell_num),:);

% real data
embryo_vid_cell = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    embryo_vid_cell{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    %embryo_vid{i} = imresize3(embryo_vid{i},round([h/5, w/5, z]));
end
embryo_vid = cat(4, embryo_vid_cell{:});
dets = cell(numel(tif_files),1);
for i=1:numel(tif_files)
    s = regionprops3(org_refine_res{i}, 'Centroid');
    dets{i} = [s.Centroid(:,2) s.Centroid(:,1) s.Centroid(:,3)];
end
movieInfo_tmp = tree2tracks_cell(org_refine_res, false);
movieInfo_tmp.tracks = cell(length(track_head), 1);
for i = 1:length(track_head)
    tr_id = track_head(i);
    if isnan(tr_id)
        movieInfo_tmp.tracks{i} = track_head;
    else
        movieInfo_tmp.tracks{i} = movieInfo.tracks{tr_id};
    end
end
% save the bad results
video_folder = 'way_2';
drawTracksStart1st(embryo_vid, dets,fullfile(res_folder,disp_folder,video_folder),...
    movieInfo_tmp);
% save ground truth
video_folder = 'way_2_gt';
movieInfo_tmp.tracks = gt_tracks;%gt_tracks(nanmax(ratio_mat,[],2)<0.8);

drawTracksStart1st(embryo_vid, dets,fullfile(res_folder,disp_folder,video_folder),...
    movieInfo_tmp);
% save one track
tr_id = 1;
save_tmp_disp(movieInfo, org_refine_res, tr_id, res_folder);
% save color cells
color_im = display_cl_frame(movieInfo, org_refine_res);
color_im = cat(5, color_im{:});
color_im = color_im*255;
drawTracksStart1st(embryo_vid, dets,fullfile(res_folder,disp_folder,...
    'color',video_folder),movieInfo_tmp, color_im);

%% use linking results to refine segmentation
movieInfo_tmp.tracks = movieInfo.tracks;
movieInfo_tmp.nei = movieInfo.nei;
tic;
refine_res_new = track2seg(movieInfo_tmp, org_refine_res, embryo_vid_cell);
toc;

% save the iterative results
movieInfo_tmp = tree2tracks_cell(refine_res_new, false);
movieInfo_tmp.tracks = cell(length(track_head), 1);
for i = 1:length(track_head)
    tr_id = track_head(i);
    if isnan(tr_id)
        movieInfo_tmp.tracks{i} = track_head;
    else
        movieInfo_tmp.tracks{i} = movieInfo.tracks{tr_id};
    end
end
dets = cell(numel(tif_files),1);
for i=1:numel(tif_files)
    s = regionprops3(refine_res_new{i}, 'Centroid');
    dets{i} = [s.Centroid(:,2) s.Centroid(:,1) s.Centroid(:,3)];
end
video_folder = 'way_2_1';
if ~exist(fullfile(res_folder,disp_folder,video_folder),'dir')
    mkdir(fullfile(res_folder,disp_folder,video_folder));
end
color_im = display_cl_frame(movieInfo, refine_res_new);
color_im = cat(5, color_im{:});
color_im = color_im*255;
drawTracksStart1st(embryo_vid, dets,fullfile(res_folder,disp_folder,video_folder),...
    movieInfo_tmp, color_im);

% file_names = cell(numel(tif_files),1);
% for i=1:numel(tif_files)
%     file_names{i} = ['Res_syn_',tif_files(i).name(1:end-4)];
% end
% display_cl_frame(movieInfo, refine_res_new, fullfile(res_folder, disp_folder, video_folder),...
%     file_names);