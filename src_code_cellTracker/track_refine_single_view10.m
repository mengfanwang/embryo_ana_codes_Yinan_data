% NOTE: script created 10022022 for processing Yinan's data with single
% view.
clc;clear;dbstop if error;
% modify term
% timepts_to_process wei_refine_res_folder save_folder 
% data_folder cell_seg_res_folder 
% g.translation_path q.saveInterMediateRes = true; 
%% define paths
cropped_flag = false;
timepts_to_process = [];
if isunix
    home_folder = fullfile('/home', getenv('USER'));
    addpath('/home/ccwang/Dropbox/cc_ImHandle');
    if cropped_flag
        data_on_insync_folder = '/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/Mengfan_data_10012022/cropped/';
        save_folder = fullfile(home_folder, data_on_insync_folder);
    else        
    
        addpath('src_code_matlab');
        timepts_to_process = generate_tps_str(0:191);
%         timepts_to_process = generate_tps_str(121:160);
        wei_refine_res_folder = "/ssd1/Mengfan/sameViewDetection_10/Wei_refine_res";
        save_folder = '/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191';
%         save_folder = '/work/Nova/embryo_res_folder/mengfan_data_res/view10_0323_000_014';
    end
    if ~exist(save_folder,'dir')
        mkdir(save_folder);
    end
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
    save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
end
addpath('dt');
addpath(genpath('CINDA/'));
addpath('src_code_matlab');
addpath('src_code_cellTracker');
addpath('src_code_cellSegment');
addpath('debug_func/');
% addpath('../../../ParticleTracking/src_code/');
% addpath('../../../ParticleTracking/src_uTrack/');
% addpath('/home/ccwang/Insync/ccwang@vt.edu/Google Drive/Projects/ParticleTracking/src_code/');
% addpath('/home/ccwang/Insync/ccwang@vt.edu/Google Drive/Projects/ParticleTracking/src_uTrack/');
addpath('/home/mengfan/ForExecute/ParticleTracking/src_code/');
addpath('/home/mengfan/ForExecute/ParticleTracking/src_uTrack/');
% set saving folder
if cropped_flag
    data_folder = fullfile(save_folder, '/imgs');
    cell_seg_res_folder = save_folder;
else

     data_folder = '/ssd1/Mengfan/sameViewFusion_10';
     cell_seg_res_folder = '/ssd1/Mengfan/sameViewDetection_10';

end
if ~exist(cell_seg_res_folder,'dir')
    disp("---> ERROR! None-exist segmentation results!");
    return;
end

%% data preparation
profile off;
% profile on;
% First downsample the detection results if any
sc_f = 2; % we resize the data to [h/sc_f, w/sc_f, z, t]

st_loc = [];
sz_crop = [];
% st_loc = [501, 1001, 101]
% sz_crop = [250, 250, 50]
gt_mat_org = {};

% read files
tif_files = dir(fullfile(data_folder, '*.tif'));
if ~isempty(timepts_to_process)
    tif_files(numel(timepts_to_process)+1:end) = [];
    for f = 1:numel(timepts_to_process)
        tif_files(f).name = timepts_to_process(f) + '.tif';
    end
end
file_num = numel(tif_files);
% load original data and detection results and ground truth
if ~exist('refine_res','var') % segmentation results
    if ~exist(fullfile(cell_seg_res_folder, 'synQuant_refine_res'), "dir")
        if ~exist(fullfile(cell_seg_res_folder, 'synQuant_refine_res_4d_v9.mat'), "file")
            load(fullfile(cell_seg_res_folder, 'synQuant_refine_res.mat'),'refine_res', 'threshold_res');
        else
            load(fullfile(cell_seg_res_folder, 'synQuant_refine_res_4d_v9.mat'),'refine_res', 'threshold_res');
        end
    else
        % the files are saved in dir "synQuant_refine_res"
        if exist('wei_refine_res_folder', 'var')
            [refine_res, threshold_res] = matfiles2cell_scaling(wei_refine_res_folder, ...
                'synQuant_refine_res', timepts_to_process, sc_f, st_loc, sz_crop);
        else
            [refine_res, threshold_res] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'synQuant_refine_res'), ...
                'synQuant_refine_res', timepts_to_process, sc_f, st_loc, sz_crop);
        end
    end
end
if file_num ~= numel(refine_res)
    warning("The number of segmentation results does not equal to number of files!");
end
% [h, w, z] = size(refine_res{1});
% org_refine_res = refine_res(1:file_num);
% org_threshold_res = threshold_res(1:file_num);
% [refine_res, ~, ~, threshold_res, ~, ~] = data_scaling(sc_f, st_loc, ...
%     sz_crop, refine_res(1:file_num), {}, {}, threshold_res(1:file_num), {}, {});

scale_term = 500;
embryo_vid = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid_temp{1} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid_temp{1} = 255*embryo_vid_temp{1}./scale_term;
    [~, embryo_vid_temp, ~, ~, ~, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, embryo_vid_temp, {}, {}, {}, {});
    embryo_vid{i} = embryo_vid_temp{1};
end
clear embryo_vid_temp


if ~exist('varMap','var') % segmentation results
    if ~exist(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), "dir")
        load(fullfile(cell_seg_res_folder, 'varianceMap.mat'),'varMap');
    else
        [varMaps, ~] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'varianceMap'), ...
            'varianceMap', timepts_to_process, sc_f, st_loc, sz_crop);
    end
end
% [~, ~, ~, ~, varMaps, ~] = data_scaling(sc_f, st_loc, ...
%     sz_crop, {}, {}, {}, {}, varMap(1:file_num), {});
% clear varMap

if ~exist('eig_res_2d','var') % segmentation results
    if ~exist(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), "dir")
        load(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res.mat'),...
            'eig_res_2d', 'eig_res_3d');
    else
%         [eig_res_2d, eig_res_3d] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), ...
%             'synQuant_priCvt_res', timepts_to_process, sc_f, st_loc, sz_crop);
        [eig_res_2d, eig_res_3d] = matfiles2cell_scaling(fullfile(cell_seg_res_folder, 'Wei_priCvt_res'), ...
            'Wei_priCvt_res', timepts_to_process, sc_f, st_loc, sz_crop);
    end
end
eigMaps = cell(file_num,1);
for i=1:file_num
    eigMaps{i} = cell(2,1);
    eigMaps{i}{1} = eig_res_2d{i};
    eigMaps{i}{2} = eig_res_3d{i};
    eig_res_2d{i} = [];
    eig_res_3d{i} = [];
end
% [~, ~, ~, ~, ~, eigMaps] = data_scaling(sc_f, st_loc, ...
%     sz_crop, {}, {}, {}, {}, {}, eigMaps);
% clear eig_res_2d eig_res_3d

%% parameter setting
g = graphPara_cell(sum(cellfun(@(x) max(x(:)), refine_res)));%1
% tform_name = "tform_050-149_11_translation.mat";
% g.translation_path = fullfile(cell_seg_res_folder, tform_name);
q = initial_q(sc_f, true);

g.timepts_to_process = timepts_to_process;
g.translation_path = '/work/Mengfan/Embryo/Registration/nonrigid_result_view10_l5_s1_linear';
% g.translation_path = '/work/Mengfan/Embryo/Registration/nonrigid_result_view11_l5_s1_linear_debug';
g.driftInfo.grid_size = 32;    % vector numbers each dim, currently cube only
g.driftInfo.batch_size = [30 30 6];
% g.resolution = [1 1 5.86];

if ~all(g.driftInfo.grid_size .*g.driftInfo.batch_size == [960 960 192])
    error('Incorrect drift info.');
end
if g.applyDrift2allCoordinate
    error('Non-rigid version doesn''t accept the option.');
end
[y_batch, x_batch, z_batch] = meshgrid(0:g.driftInfo.grid_size+1);
if ~isempty(st_loc)  % adjust if crop
    y_batch = y_batch*g.driftInfo.batch_size(2) + 0.5 - g.driftInfo.batch_size(2)/2 - st_loc(2)/sc_f;
    x_batch = x_batch*g.driftInfo.batch_size(1) + 0.5 - g.driftInfo.batch_size(1)/2 - st_loc(1)/sc_f;
    z_batch = z_batch*g.driftInfo.batch_size(3) + 0.5 - g.driftInfo.batch_size(3)/2 - st_loc(3);
else
    y_batch = y_batch*g.driftInfo.batch_size(2) + 0.5 - g.driftInfo.batch_size(2)/2;
    x_batch = x_batch*g.driftInfo.batch_size(1) + 0.5 - g.driftInfo.batch_size(1)/2;
    z_batch = z_batch*g.driftInfo.batch_size(3) + 0.5 - g.driftInfo.batch_size(3)/2;
end
g.driftInfo.y_batch = y_batch;
g.driftInfo.x_batch = x_batch;
g.driftInfo.z_batch = z_batch;  % vector locations in original image with one padding

%% start tracking
q.saveInterMediateRes = true;      % get all resuts
if length(timepts_to_process) > 30
    q.saveInterMediateRes = false;
end
q.save_folder = save_folder;
diary(fullfile(save_folder, 'log'));

[movieInfo, movieInfoAll, out_refine_res, refine_resAll,...
    threshold_res, threshold_resAll] = .... 
    mcfTracking_cell(refine_res, embryo_vid, threshold_res, ...
    varMaps, eigMaps, g, q);
% profile viewer;

%% save results
%display_rgb_time_course(out_refine_res);

if q.saveInterMediateRes
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo',...
        'movieInfoAll','-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res',...
        'refine_resAll','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), 'threshold_res',...
        'threshold_resAll','-v7.3');
else
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo', '-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), 'threshold_res','-v7.3');
end

%% redefine result folder as a tmp folder
%save_folder = fullfile(home_folder, '/Desktop/embryo_res_folder/mengfan_data_res/Mengfan_data_10012022/');
%% save results to mastodon
mastodon_dir = fullfile(save_folder, 'mastodon');
if ~exist(mastodon_dir)
    mkdir(mastodon_dir);
end
addpath('TGMM_wrapper/');
mat2tgmm(movieInfo, fullfile(mastodon_dir, 'tgmm_format'));
tif2bdv(data_folder, fullfile(mastodon_dir, 'embryo_data_h5'), timepts_to_process, st_loc, sz_crop);
%% stop here
return;
%% save seg data for Yinan
seg_res_dir = fullfile(save_folder, 'seg_before_refine');
if ~exist(track_res_dir)
    mkdir(track_res_dir);
end
for i = 1:numel(tif_files)  
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    save_seg_masks(org_im/2000, refine_res{i}, fullfile(seg_res_dir, ...
        tif_files(i).name(1:length(tif_files(i).name)-4)))
end

%% save track data for Yinan
track_res_dir = fullfile(save_folder, 'track_res');
if ~exist(track_res_dir)
    mkdir(track_res_dir);
end
l = cellfun(@length, movieInfo.tracks);
[track_length, idx] = sort(l, 'descend');
save_cnt = 1;
for i = 1:numel(movieInfo.tracks)
    track_id = idx(i);
    if movieInfo.frames(movieInfo.tracks{track_id}(1)) == 1
        saved = save_maxProj_track(movieInfo, track_id, out_refine_res,...
            embryo_vid, track_res_dir, save_cnt);
        if saved
            save_cnt = save_cnt + 1;
        end
    end
end