% NOTE: script created 10022022 for processing Yinan's data with single
% view.
%% define paths
cropped_flag = false;
timepts_to_process = [];
if isunix
    home_folder = fullfile('/home', getenv('USER'));
    addpath(fullfile(home_folder, 'Dropbox/cc_ImHandle/'));
    if cropped_flag
        data_on_insync_folder = '/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/Mengfan_data_10012022/cropped/';
        save_folder = fullfile(home_folder, data_on_insync_folder);
    else
        % mengfan's 240-249 frame whole data
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_10012022');
        % mengfan's 230-239 frame whole data
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_230_239_10072022');
        % mengfan's 230-249 frame whole data
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_230_249_10072022');
        % 50-69 frames
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_50_69_10072022');
        % 66-88 frames: NOTE: other than |refine_res| and |threshold_res|,
        % the other structures contain 66-89 time points. This is due to
        % mengfan's code problem.
        % save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_66_88_10072022');
        
        % 50-94 frames
        %timepts_to_process = generate_tps_str(50:94);
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_50_94_1022');
        
        % 90-129 frames
        %timepts_to_process = generate_tps_str(90:129);
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_90_129_1022');
        
        % 50-149 frames
        %timepts_to_process = generate_tps_str(50:149);
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_50_149_1023');
        
        % 72-91 frames
        %st_t = 72;
        %end_t = 91;
        %timepts_to_process = generate_tps_str(st_t:end_t);
        %save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder/mengfan_data_res/mengfan_single_view_72_91_1116');
        
        % 50-149 frames using mengfan's registration
        timepts_to_process = generate_tps_str(50:149);
        wei_refine_res_folder = "/work/public/sameViewFusion/sameViewDetection_050-149_11_v6/Wei_refine_res";
        save_folder = fullfile(home_folder, ...
            'Desktop/embryo_res_folder/mengfan_data_res/result_50_149_1203_validation2');
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
addpath('../../../ParticleTracking/src_code/');
addpath('../../../ParticleTracking/src_uTrack/');
% set saving folder
if cropped_flag
    data_folder = fullfile(save_folder, '/imgs');
    cell_seg_res_folder = save_folder;
else
    % mengfan's 240-249 frame whole data
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_240-249_08';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_240-249_08';
    % mengfan's 230-249 frame whole data
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_230-239_08';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_230-239_08';
    % mengfan's 230-249 frame whole data
%     data_folder = '/work/public/sameViewFusion/TwentyFrames/sameViewFusion_230-249_08';
%     cell_seg_res_folder = '/work/public/sameViewFusion/TwentyFrames/sameViewDetection_230-249_08';

    % 50-69 frame whole data
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_050-069_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-069_11';
    % 66-88 frame whole data
%     data_folder = '/work/public/sameViewFusion/merged_batch/sameViewFusion_066-088_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/merged_batch/sameViewDetection_066-088_11';
    
    % 50-94 frames
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_050-129_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-129_11';
    % 90-129 frames --> same as 50-94
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_090-129_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11';
    
%     % 50-149 frames 
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_050-149_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11';

    % 72-91 frames --> same as 50-149
%     data_folder = '/work/public/sameViewFusion/sameViewFusion_050-149_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11';
    
    % 50-149 frames using Wei's new results
     data_folder = '/work/public/sameViewFusion/sameViewFusion_050-149_11';
     cell_seg_res_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11';

%----------------------------------TODO----------------------------------------%
    % 50-74 frame whole data
%     data_folder = '/work/public/sameViewFusion/merged_batch/sameViewFusion_050-74_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/merged_batch/sameViewDetection_050-74_11';
    % 70-88 frame whole data
%     data_folder = '/work/public/sameViewFusion/merged_batch/sameViewFusion_070-89_11';
%     cell_seg_res_folder = '/work/public/sameViewFusion/merged_batch/sameViewDetection_070-89_11';
end
if ~exist(cell_seg_res_folder,'dir')
    disp("---> ERROR! None-exist segmentation results!");
    return;
end

%% data preparation
profile off;
profile on;
% First downsample the detection results if any
sc_f = 2; % we resize the data to [h/sc_f, w/sc_f, z, t]

st_loc = [];
sz_crop = [];
% st_loc = [251, 1, 1];
% sz_crop = [200, 250, z];
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
            [refine_res, threshold_res] = matfiles2cell(wei_refine_res_folder, ...
                'synQuant_refine_res', timepts_to_process);
        else
            [refine_res, threshold_res] = matfiles2cell(fullfile(cell_seg_res_folder, 'synQuant_refine_res'), ...
                'synQuant_refine_res', timepts_to_process);
        end
    end
end
if file_num ~= numel(refine_res)
    warning("The number of segmentation results does not equal to number of files!");
end
[h, w, z] = size(refine_res{1});
% org_refine_res = refine_res(1:file_num);
% org_threshold_res = threshold_res(1:file_num);
[refine_res, ~, ~, threshold_res, ~, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, refine_res(1:file_num), {}, {}, threshold_res(1:file_num), {}, {});

scale_term = 1000;
embryo_vid = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid{i} = 255*embryo_vid{i}./scale_term;
end
[~, embryo_vid, ~, ~, ~, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, embryo_vid, {}, {}, {}, {});

if ~exist('varMap','var') % segmentation results
    if ~exist(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), "dir")
        load(fullfile(cell_seg_res_folder, 'varianceMap.mat'),'varMap');
    else
        [varMap, ~] = matfiles2cell(fullfile(cell_seg_res_folder, 'varianceMap'), ...
            'varianceMap', timepts_to_process);
    end
end
[~, ~, ~, ~, varMaps, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, {}, {}, {}, varMap(1:file_num), {});
clear varMap

if ~exist('eig_res_2d','var') % segmentation results
    if ~exist(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), "dir")
        load(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res.mat'),...
            'eig_res_2d', 'eig_res_3d');
    else
        [eig_res_2d, eig_res_3d] = matfiles2cell(fullfile(cell_seg_res_folder, 'synQuant_priCvt_res'), ...
            'synQuant_priCvt_res', timepts_to_process);
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
[~, ~, ~, ~, ~, eigMaps] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, {}, {}, {}, {}, eigMaps);
% clear eig_res_2d eig_res_3d

%% parameter setting
g = graphPara_cell(sum(cellfun(@(x) max(x(:)), refine_res)));%1
tform_name = "tform_050-149_11_translation.mat";
g.translation_path = fullfile(cell_seg_res_folder, tform_name);
q = initial_q(sc_f, true);

%% start tracking
[movieInfo, movieInfoAll, out_refine_res, refine_resAll,...
    threshold_res, threshold_resAll] = ...
    mcfTracking_cell(refine_res, embryo_vid, threshold_res, ...
    varMaps, eigMaps, g, q);
profile viewer;

%% save results
%display_rgb_time_course(out_refine_res);

if q.saveInterMediateRes
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo',...
        'movieInfoAll','-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res',...
        'refine_resAll','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), threshold_res,...
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
tif2bdv(data_folder, fullfile(mastodon_dir, 'embryo_data_h5'));
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