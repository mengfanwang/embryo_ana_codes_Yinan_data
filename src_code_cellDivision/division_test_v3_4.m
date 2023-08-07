clc;clear;close all;
% v3.3 = v3.1 + v3.2 two-step module
% first: check two component detections. If size ratio < 3, it's a division.
% second: use v3.1 pipline to detect remaining cases
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0504_000_191/movieInfo.mat');
%% get ground truth
addpath ../gt_measure/
path = '/work3/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(path, 1);
% link table: s t
% spot table: id t y x z
link_table = link_table(:,1:2);
spot_table = spot_table(:,1:5);

% convert to tracking result corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1/2 1/2 1/5.86] + 1;

%% find all detections and candidate divisions in movieInfo
t = 192;
data_size = [960 960 181];
cell_num = size(spot_table,1);
detection_list = nan(cell_num,1);
div_candiList = cell(cell_num,1);
fprintf('Find cells and candidate divisions in movieInfo...\n'); 
% get heads of tracks
head_list = nan(length(movieInfo.tracks),1);
for ii = 1:length(movieInfo.tracks)
    if length(movieInfo.tracks{ii}) >= 10
        head_list(ii) = movieInfo.tracks{ii}(1);
    end
end
Ci = abs(movieInfo.Ci(1));
for ii = 1:cell_num
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, cell_num);
    end
    cell_loc = round(spot_table(ii, 3:5));
    cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
    cell_time = spot_table(ii, 2);
    ind_candi = find(movieInfo.frames == cell_time);
    ind_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(ind_candi)));
    if isempty(ind_match)
        warning('Missing cell.');
    else
        detection_id = ind_match + sum(movieInfo.n_perframe(1:cell_time-1));
        detection_list(ii) = detection_id;
        % if have child already, it should be small
        child_small_flag = true;
        child_id = movieInfo.kids{detection_id};
        if ~isempty(child_id)
            if length(movieInfo.voxIdx{child_id}) >= ...
                length(movieInfo.voxIdx{detection_id}) * 0.8   
                child_small_flag = false;
            end
        end
        if cell_time ~= t && child_small_flag
            % find the 5-nearest neighbor
            div_candi = findNeighbor(movieInfo, detection_id);
            % only heads are included   
            div_candi = div_candi(ismember(div_candi,head_list));
            if isempty(div_candi)
                continue;
            end

            div_candiList{ii} = div_candi;
        end
    end
end
                                            

%% convert candiList to division cases
% should statisify three creteria
% 1.every child can happen only once
% 2.parent must be the nearest one to the center of the pair of children
% 3.children must split
div_num = 0;
detect_divPair = zeros(10000,3);
for ii = 1:length(detection_list)
    parent_id = detection_list(ii);
    if isnan(parent_id) || isempty(div_candiList{ii})
        continue;
    end
    if ~isempty(movieInfo.kids{parent_id})     % case1: one child exists
        child1_id = movieInfo.kids{parent_id};
        child_time = movieInfo.frames(child1_id);
        gap_flag = false(length(div_candiList{ii}),1);
        dist2parent = inf(length(div_candiList{ii}),1);
        for jj = 1:length(div_candiList{ii})
            child2_id = div_candiList{ii}(jj);
            gap_flag(jj) = division_gap_test(movieInfo, zeros(data_size), child1_id, child2_id);
            frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
                movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
%             dist2parent(jj) = ovDistanceRegion(movieInfo.vox{parent_id},...
%                 [movieInfo.vox{child1_id}; movieInfo.vox{child2_id}], frame_shift);
            dist2parent(jj) = norm(mean(movieInfo.vox{parent_id}) + frame_shift -...
                mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}]));
        end
        [dist2parent_min, minIdx] = min(dist2parent);
        if gap_flag(minIdx)
            div_num = div_num + 1;
            detect_divPair(div_num,:) = [parent_id child1_id div_candiList{ii}(minIdx)];
        end
        
    elseif length(div_candiList{ii}) >= 2
        children_pair = nchoosek(div_candiList{ii}, 2);
        child_time = movieInfo.frames(parent_id)+1;
        gap_flag = false(size(children_pair,1),1);
        dist2parent = inf(size(children_pair,1),1);
        for jj = 1:size(children_pair,1)
            child1_id = children_pair(jj,1);
            child2_id = children_pair(jj,2);
            gap_flag(jj) = division_gap_test(movieInfo, zeros(data_size), child1_id, child2_id);
            frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
                movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
%             dist2parent(jj) = ovDistanceRegion(movieInfo.vox{parent_id},...
%                 [movieInfo.vox{child1_id}; movieInfo.vox{child2_id}], frame_shift);
            dist2parent(jj) = norm(mean(movieInfo.vox{parent_id}) + frame_shift -...
                mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}]));
        end
        [dist2parent_min, minIdx] = min(dist2parent);
        if gap_flag(minIdx)
            div_num = div_num + 1;
            detect_divPair(div_num,:) = [parent_id children_pair(minIdx,1) children_pair(minIdx,2)];
        end
    end
end
detect_divPair = detect_divPair(1:div_num,:);
detect_divPair = unique(detect_divPair, 'rows');

%% for each pair of the children, mapping their center to the nearest parent
multi_child_list = tabulate(reshape(detect_divPair(:,2:3), [] ,1));
multi_child_list = find(multi_child_list(:,2) > 1);
parent_remove_flag = zeros(size(detect_divPair,1),1);
for ii = 1:size(multi_child_list,1)
    parent_loc = find(any(ismember(detect_divPair(:,2:3), multi_child_list(ii)),2));
    dist2parent = inf(size(parent_loc));
    child_time = movieInfo.frames(multi_child_list(ii));
    for jj = 1:length(parent_loc)
        parent_id = detect_divPair(parent_loc(jj),1);
        child1_id = detect_divPair(parent_loc(jj),2);
        child2_id = detect_divPair(parent_loc(jj),3);
        frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
            movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
%         dist2parent(jj) = ovDistanceRegion(movieInfo.vox{parent_id},...
%             [movieInfo.vox{child1_id}; movieInfo.vox{child2_id}], frame_shift);
        dist2parent(jj) = norm(mean(movieInfo.vox{parent_id}) + frame_shift -...
            mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}]));
    end
     [~, minIdx] = min(dist2parent);
    parent_remove_flag(parent_loc) = 1;
    parent_remove_flag(parent_loc(minIdx)) = 0;    

end
detect_divPair(logical(parent_remove_flag), :) = [];

%% merge 3.1 results and 3.2 results
load('detect_splitPair.mat');
% remove divPair cases happened in splitPair
detect_remove_flag = true(size(detect_divPair,1),1);
for ii = 1:size(detect_divPair,1)
    if any(ismember(detect_divPair(ii,:), detect_splitPair(:)))
        detect_remove_flag(ii) = false;
    end
end
detect_allPair = [detect_divPair(detect_remove_flag,:); detect_splitPair;];

%% check division detection accuracy
[TP, TN, FP, gt_divPair, gt_code_list, detect_code_list] = ...
    div_acc(detect_allPair, G, movieInfo, spot_table, data_size);
fprintf('TP: %d  TN: %d  FP: %d\n', TP, TN, FP);
% 5 47 69 75 77 81: wrong association
% 10 24 68 83 87: sudden appear cell
% 52: neighbor borken 
% 111: FP deteciton
% 112: over-merge

%% FP investigation
gt_fnPair = gt_divPair(gt_code_list~=1,:);
cell_match_flag = nan(size(gt_fnPair,1), 3);
size_ratio = nan(size(gt_fnPair,1), 2);
track_loc = nan(size(gt_fnPair,1), 2);
track_length = nan(size(gt_fnPair,1), 2);
parent_track_flag = false(size(gt_fnPair,1), 2);
not_neighbor_flag = false(size(gt_fnPair,1), 2);

not_nearest_parent_flag = false(size(gt_fnPair,1), 2);
dist2parent = inf(size(gt_fnPair,1), 2);
dist2nearest = inf(size(gt_fnPair,1), 2);
for ii = 1:size(gt_fnPair,1)
    parent_id = findMatchCell(movieInfo, gt_fnPair(ii,1), gt_fnPair(ii,4));
    child1_id = findMatchCell(movieInfo, gt_fnPair(ii,2), gt_fnPair(ii,4)+1);
    child2_id = findMatchCell(movieInfo, gt_fnPair(ii,3), gt_fnPair(ii,4)+1);
    cell_match_flag(ii,:) = [parent_id child1_id child2_id];
    if ~isnan(parent_id)
        parent_track = movieInfo.particle2track(parent_id,1);
        neighbors = findNeighbor(movieInfo, parent_id);
        for jj = 1:2
            if jj == 1
                child_id = child1_id;
            elseif jj == 2
                child_id = child2_id;
            end
            if ~isnan(child_id)
                size_ratio(ii,jj) = length(movieInfo.voxIdx{child_id})...
                    ./ length(movieInfo.voxIdx{parent_id});
                child_track = movieInfo.particle2track(child_id,1);
                not_neighbor_flag(ii,jj) = ~ismember(child_id, neighbors);
                if child_track == parent_track
                    parent_track_flag(ii,jj) = 1;
                end
                if isnan(child_track)            % single detection
                    track_loc(ii,jj) = 1;
                    track_length(ii,jj) = 1;
                else
                    track_loc(ii,jj) = find(movieInfo.tracks{child_track}==child_id);
                    track_length(ii,jj) = length(movieInfo.tracks{child_track}) - track_loc(ii,jj) + 1;
                end

                % the parent should be the nearest as well
                parent_candi = findParent(movieInfo, child_id);      
                dist2parent = inf(size(parent_candi));
                frame_shift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(child_id,:),...
                    gt_fnPair(ii,4), gt_fnPair(ii,4)+1, movieInfo.drift, movieInfo.driftInfo);
                for kk = 1:length(parent_candi)
                    dist2parent(kk) = ovDistanceRegion(movieInfo.vox{parent_candi(kk)},...
                    movieInfo.vox{child_id}, frame_shift);
                end
                [dist2parent_min, minIdx] = min(dist2parent);
                dist2nearest(ii,jj) = dist2parent_min;
                if ismember(parent_id, parent_candi)
                    dist2parent(ii,jj) = dist2parent(parent_candi == parent_id);
                end
                if parent_candi(minIdx) ~= parent_id || dist2parent_min >= 33
                    not_nearest_parent_flag(ii,jj) = 1;
                end


            end
        end
    end
end

fprintf('Missing cell #:%d %f\n', sum(any(isnan(cell_match_flag)')), ...
    sum(any(isnan(cell_match_flag)'))/size(gt_fnPair,1));

too_large_flag = size_ratio>0.8;
fprintf('Large child #:%d %f\n', sum(any(too_large_flag')), ...
    sum(any(too_large_flag'))/size(gt_fnPair,1));

over_merge_flag = cell_match_flag(:,2) == cell_match_flag(:,3);
for ii = 1:size(gt_fnPair,1)
    if over_merge_flag(ii)
        if ismember(cell_match_flag(ii,:), detect_splitPair, 'rows')
            over_merge_flag(ii) = 0;
        end
    end
end
fprintf('Over merge #:%d %f\n', sum(over_merge_flag), ...
    sum(over_merge_flag)/size(gt_fnPair,1));

not_head_flag = track_loc;
not_head_flag(parent_track_flag) = nan;
not_head_flag = not_head_flag > 1;
fprintf('Not head#:%d %f\n', sum(any(not_head_flag,2)), ...
    sum(any(not_head_flag,2))/size(gt_fnPair,1));

too_short_flag = track_length;
too_short_flag(parent_track_flag) = nan;
too_short_flag = too_short_flag < 10;
fprintf('Too short#:%d %f\n', sum(any(too_short_flag,2)), ...
    sum(any(too_short_flag,2))/size(gt_fnPair,1));

fprintf('Not neighbor#:%d %f\n', sum(any(not_neighbor_flag,2)), ...
    sum(any(not_neighbor_flag,2))/size(gt_fnPair,1));

fprintf('Not nearest parent#:%d %f\n', sum(any(not_nearest_parent_flag,2)), ...
    sum(any(not_nearest_parent_flag,2))/size(gt_fnPair,1));

sum(over_merge_flag | any(too_large_flag, 2) | any(isnan(cell_match_flag),2) ...
    | any(not_head_flag,2) | any(too_short_flag,2) | any(not_neighbor_flag,2) | any(not_nearest_parent_flag,2) ) 
%% convert FP to movieInfo format for visiualization
fp_pair = detect_allPair(detect_code_list~=1,:);
fp_tracks = cell(size(fp_pair*2,1), 1);
movieInfo_temp = movieInfo;
movieInfo_temp.parents = cell(size(movieInfo.parents));
movieInfo_temp.kids = cell(size(movieInfo.kids));
for ii = 1:size(fp_pair,1)
    fp_tracks{ii*2-1} = fp_pair(ii,[1 2]);
    movieInfo_temp.parents{fp_pair(ii,2)} = fp_pair(ii,1);
    movieInfo_temp.kids{fp_pair(ii,1)} = fp_pair(ii,2);
    
    fp_tracks{ii*2} = fp_pair(ii,[1 3]);
    movieInfo_temp.parents{fp_pair(ii,3)} = fp_pair(ii,1);
    movieInfo_temp.kids{fp_pair(ii,1)} = fp_pair(ii,3);
end
movieInfo_temp.tracks = fp_tracks;
mat2tgmm(movieInfo_temp, '/work/Nova/embryo_res_folder/mengfan_data_res/temp');

%%
function gap_flag = division_gap_test(movieInfo, embryo_vid, child1, child2)
    % get local region
    [vidOut, vIdx, ~, ~] = crop3D(embryo_vid, [movieInfo.voxIdx{child1}; movieInfo.voxIdx{child2}], [0 0 0]);
    vCate = zeros(size(vidOut));
    vCate(ismember(vIdx, movieInfo.voxIdx{child1})) = 1;
    vCate(ismember(vIdx, movieInfo.voxIdx{child2})) = 2;
%     % get convex hull
%     [x_ind, y_ind, z_ind] = ind2sub_direct(size(vCate), find(vCate));
%     [x_test, y_test, z_test] = ind2sub_direct(size(vCate), find(ones(size(vCate))));
%     hull_flag = inhull([x_test, y_test, z_test],[x_ind, y_ind, z_ind],...
%          convhulln([x_ind y_ind z_ind]));
%     vCate(hull_flag) = 3;
%     vCate(ismember(vIdx, movieInfo.voxIdx{child1})) = 1;
%     vCate(ismember(vIdx, movieInfo.voxIdx{child2})) = 2;
%     % z_test
%     fg_vals = vidOut(vCate==1 | vCate==2);
%     bg_vals = vidOut(vCate==3);
%     z_score = (mean(fg_vals) - mean(bg_vals))/std(bg_vals)*sqrt(length(fg_vals));
    a = regionprops3(vCate>0, "VoxelIdxList");
    if size(a,1) == 2
        gap_flag = true;
    else
        gap_flag = false;
    end
end

function [TP, TN, FP, gt_divPair, gt_code_list, detect_code_list] = ...
    div_acc(detect_divPair, G, movieInfo, spot_table, data_size) 
    % get ground truth pairs
    d_out = outdegree(G);
    gt_list = G.Nodes.Name(find(d_out == 2));
    gt_divPair = nan(size(gt_list,1), 4);  %[parent_idx child_idx*2 parent_time]
    for ii = 1:length(gt_list)
        gt_pair = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_list{ii}),2);
        gt_pair = {gt_list{ii}; gt_pair{1}; gt_pair{2}};
        for jj = 1:3
            cell_ind = str2num(gt_pair{jj});
            cell_ind = find(spot_table(:,1) == cell_ind);
            cell_loc = round(spot_table(cell_ind, 3:5));
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            gt_divPair(ii,jj) = cell_idx;
            if jj == 1
                gt_divPair(ii, 4) = spot_table(cell_ind, 2);
            end
        end
    end
    % add missed annotation
    gt_addPair = [53044942 51210387 55802059 40;
                  116506827 113737225 117436109 60;
                  89013364 89007599 86248573 80;
                  66755697 68589298 64919216 125;
                  82375979 83288939 83303339 5;
                  78923752 80763122 79851106 104];
    gt_divPair = [gt_divPair; gt_addPair];
    % remove one incorrect annotation
    error_loc = find(cellfun(@(x)strcmp(x, '34390'), gt_list));
    gt_list(error_loc) = [];
    gt_divPair(error_loc,:) = [];
    gt_correctPair = [41064118 41989558 41059316 17;];
    gt_divPair = [gt_divPair; gt_correctPair];
  

    gt_code_list = zeros(size(gt_divPair,1), 1);
    detect_code_list = zeros(size(detect_divPair,1),1);
    % error code:
    % 0: missing parent
    % 1: correct
    % 2: child1 incorrect
    % 3: child2 incorrect
    % 4: both incorrect 
    for ii = 1:size(gt_divPair,1)
        parent_id = findMatchCell(movieInfo, gt_divPair(ii,1), gt_divPair(ii,4));
        detect_loc = find(parent_id == detect_divPair(:,1));
        if isnan(parent_id) || isempty(detect_loc)                          % case 0: missing parent
            error_code = 0;
        else
            child1_id = findMatchCell(movieInfo, gt_divPair(ii,2), gt_divPair(ii,4)+1);
            child2_id = findMatchCell(movieInfo, gt_divPair(ii,3), gt_divPair(ii,4)+1);
            if ismember(child1_id, detect_divPair(detect_loc,2:3)) && ...
               ismember(child2_id, detect_divPair(detect_loc,2:3))          % case 1: correct
                gt_code_list(ii) = 1;
                detect_code_list(detect_loc) = 1;
                continue;
            elseif ~ismember(child1_id, detect_divPair(detect_loc,2:3)) 
                error_code = 2;                                             % case 2: child1 incorrect
            elseif ~ismember(child2_id, detect_divPair(detect_loc,2:3)) 
                error_code = 3;                                             % case 3: child2 incorrect   
            else
                error_code = 4;                                             % case 4: both incorrect
            end
        end
        gt_code_list(ii) = error_code;
        if ~isempty(detect_loc)
            detect_code_list(detect_loc) = error_code;
        end
    end

    
    TP = sum(gt_code_list == 1);
    TN = length(gt_code_list) - TP;
    FP = size(detect_divPair,1) - TP;
end

function cell_match = findMatchCell(movieInfo, cell_idx, cell_time)
    cell_candi = find(movieInfo.frames == cell_time);
    cell_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(cell_candi)));
    if isempty(cell_match)
        cell_match = nan;
    else
        cell_match = cell_match + sum(movieInfo.n_perframe(1:cell_time-1));
    end
end


function neighbors = findNeighbor(movieInfo, cur_i)
    % find the nearest neighbors in the next frame with given conditions

    max_dist = 50;
    max_nei = 5;
    resolution = [1 1 5.86];
    
    tt = movieInfo.frames(cur_i);
    curCentroid = movieInfo.orgCoord(cur_i,:);
    nextCentroid = movieInfo.orgCoord(movieInfo.frames==(tt+1),:);
    
    % candidate neighbors roughly selection
    drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt, tt+1, movieInfo.drift, movieInfo.driftInfo);
    dist_pair = pdist2((curCentroid + drift).*resolution, nextCentroid.*resolution); % not accurate but enough
    [neighbor_dist, neighbor_candidate] = sort(dist_pair);
    neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
    neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
    neighbors = neighbor_candidate(neighbor_dist < max_dist);
    neighbors = neighbors + sum(movieInfo.n_perframe(1:tt));
end

function neighbors = findParent(movieInfo, cur_i)
    % find the nearest neighbors in the next frame with given conditions

    max_dist = 50;
    max_nei = 5;
    resolution = [1 1 5.86];
    
    tt = movieInfo.frames(cur_i);
    curCentroid = movieInfo.orgCoord(cur_i,:);
    prevCentroid = movieInfo.orgCoord(movieInfo.frames==(tt-1),:);
    
    % candidate neighbors roughly selection
    drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt-1, tt, movieInfo.drift, movieInfo.driftInfo);
    dist_pair = pdist2((curCentroid + drift).*resolution, prevCentroid.*resolution); % not accurate but enough
    [neighbor_dist, neighbor_candidate] = sort(dist_pair);
    neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
    neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
    neighbors = neighbor_candidate(neighbor_dist < max_dist);
    neighbors = neighbors + sum(movieInfo.n_perframe(1:tt-2));
end
