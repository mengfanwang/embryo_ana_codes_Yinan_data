clearvars -except movieInfo embryo_vid
dbstop if error
% load data
% load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');

%% INPUT:
% movieInfo ------------------ result from tracking framework
% embryo_vid  ---------------- imaging data (may scaled before tracking)
% OUTPUT:
% detect_divPair ------------- n*4 matrix, [parent_id child_id*2]
% mengfanw@vt.edu, 08/25/2023

%% parameter setting
paras.im_resolution = [1 1 5.86/2];
paras.max_dist = 50/2;    % max cell moving distance (normal to resolution)

paras.max_nei = 5;        % max neighbor can find
paras.size_ratio = 3;     % max size ratio of two children
paras.edge_ratio = 5;     % max edge ratio of two children
paras.length_thre = 10;   % minimum track length to find a division
paras.child2parent_ratio = 0.8; % maximum child to parent ratio

paras.data_size = size(embryo_vid{1});

%% main functions
% step 1: get a part of divisions by checking components
detect_splitPair = component_divDetection(movieInfo, paras);
% step 2: get remaing divisions by rule-based method
detect_divPair = rule_based_divDetection(movieInfo, embryo_vid, paras);
detect_divPair = detect_divPair(:,1:3);   % 4 is for debug
detect_divPair = [detect_divPair; detect_splitPair];

%% build new tracks and save to mastodon file
movieInfo_temp = movieInfo;
for ii = 1:size(detect_divPair,1)
    if mod(ii,1000) == 0
        fprintf('%d/%d\n', ii, size(detect_divPair,1));
    end
    parent_id = detect_divPair(ii,1);
    child1_id = detect_divPair(ii,2);
    child2_id = detect_divPair(ii,3);

    parent_track = movieInfo.particle2track(parent_id,1);
    child1_track = movieInfo.particle2track(child1_id,1);
    child2_track = movieInfo.particle2track(child2_id,1);
    
    if parent_track == child1_track
        % parent connects to child1, now add to child2
        movieInfo_temp = divFun.mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child2_id, child2_track);
    else 
        % parents are not connected to any child
        movieInfo_temp = divFun.mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child1_id, child1_track);
        movieInfo_temp = divFun.mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child2_id, child2_track);
    end
end

%% save to mastodon file
mat2tgmm(movieInfo_temp, '/work/Nova/embryo_res_folder/mengfan_data_res/merge_division/tgmm_format');
