clc;clear;close all;
% change the evaluation: check division cases rather than edges
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0504_000_191/movieInfo.mat');
%% get ground truth
addpath ../gt_measure/
path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(path, 1);
% link table: s t
% spot table: id t y x z
link_table = link_table(:,1:2);
spot_table = spot_table(:,1:5);

% convert to tracking result corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1/2 1/2 1/5.86] + 1;

%% convert to tracks
d_in = indegree(G);
d_out = outdegree(G);

component_id_list = conncomp(G, 'Type', 'weak')';
tracks = cell(max(component_id_list),1);
link_flag = zeros(size(link_table,1),1);

for component_id = 1:max(component_id_list)
    component_id
    tracks{component_id} = {};
    track_id = 0;
    cell_ind = find(component_id_list == component_id & d_in == 0);
    s_node = G.Nodes.Name{cell_ind};
    t_ind_list = find(component_id_list == component_id & d_out == 0);
    t_node_list = G.Nodes.Name(t_ind_list);
    for ii = 1:length(t_node_list)
        track_id = track_id + 1;
        tracks{component_id}{track_id} = [];
        path = allpaths(G, s_node, t_node_list{ii});
        if length(path) > 2
            error('Multi paths!');
        end
        path = path{1};
        for jj = 1:length(path) - 1
            parent = str2num(path{jj});
            child = str2num(path{jj+1});
            link_ind = find(ismember(link_table, [parent child], 'rows'));
            if link_flag(link_ind) == 0
                if isempty(tracks{component_id}{track_id})
                    tracks{component_id}{track_id} = [parent child];
                else
                    tracks{component_id}{track_id}(end+1) = child;
                end
                link_flag(link_ind) = 1;
            end
        end
    end
end

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
            div_candi = findNeighbor(movieInfo, detection_id)...
                + sum(movieInfo.n_perframe(1:cell_time));
            % only heads are included   
            div_candi = div_candi(ismember(div_candi,head_list));
            if isempty(div_candi)
                continue;
            end

            % double check the parent should be the nearest as well
            parent_nearest_flag = zeros(size(div_candi));
            for jj = 1:length(div_candi)
                parent_candi = findParent(movieInfo, div_candi(jj))...
                    + sum(movieInfo.n_perframe(1:cell_time-1));      
                parent_cost = zeros(size(parent_candi));
                frame_shift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(div_candi(jj),:),...
                    cell_time, cell_time+1, movieInfo.drift, movieInfo.driftInfo);
                for kk = 1:length(parent_candi)
                    parent_cost(kk) = ovDistanceRegion(movieInfo.vox{parent_candi(kk)},...
                    movieInfo.vox{div_candi(jj)}, frame_shift);
                end
                [parent_cost_min, parent_min] = min(parent_cost);

                if parent_candi(parent_min) == detection_id && parent_cost_min < 33
                    parent_nearest_flag(jj) =1;
                end
            end
            div_candi = div_candi(logical(parent_nearest_flag));
            if isempty(div_candi)
                continue;
            end

            div_candiList{ii} = div_candi;
        end
    end
end

%% gap test between two children -- loose creteria
% if one child exists, test it with other children
% elseif only one child, keep it
% elseif multiple children, test all pairs and remove if all fails
div_candiList_old = div_candiList;
for ii = 1:length(detection_list)
    parent_id = detection_list(ii);
    if isnan(parent_id) || isempty(div_candiList_old{ii})
        continue;
    end
    if ~isempty(movieInfo.kids{parent_id})     % case1: one child exists
        child1 = movieInfo.kids{parent_id};
        child_time = movieInfo.frames(child1);
        for jj = 1:length(div_candiList_old{ii})
            child2 = div_candiList_old{ii}(jj);
            z_score = division_gap_test(movieInfo, embryo_vid{child_time}, child1, child2);
%             z_threshold = -norminv(0.01/length(movieInfo.vox)); % FDR control
            if z_score == 0
                div_candiList{ii}(div_candiList{ii} == child2) = [];
            end
        end
    end
end


%% convert candiList to division cases
div_num = 0;
detect_divPair = zeros(1000,3);
for ii = 1:length(detection_list)
    parent_id = detection_list(ii);
    if isnan(parent_id) || isempty(div_candiList{ii})
        continue;
    end
    if ~isempty(movieInfo.kids{parent_id})     % case1: one child exists
        if length(div_candiList{ii}) > 1
            error('Multi candidates!');
        end
        child1 = movieInfo.kids{parent_id};
        child2 = div_candiList{ii};
    else
        if length(div_candiList{ii}) == 2
            child1 = div_candiList{ii}(1);
            child2 = div_candiList{ii}(2);
        elseif length(div_candiList{ii}) > 2
            error('Multi candidates!');
        end
    end
    div_num = div_num + 1;
    detect_divPair(div_num,:) = [parent_id child1 child2];
end
detect_divPair = detect_divPair(1:div_num,:);

% for each pair of the parents, mapping their center to the nearest parent
child_list = unique(detect_divPair(:,2:3), 'rows');
[~,child_fre] = mode(child_list,'all');
if child_fre > 1
    error('Multiple kids kept!');
end
parent_remove_flag = zeros(size(detect_divPair,1),1);
for ii = 1:size(child_list,1)
    parent_loc = find(ismember(detect_divPair(:,2:3), child_list(ii,:), 'rows'));
    if length(parent_loc) > 1
        parent_cost = zeros(size(parent_loc));
        child_time = movieInfo.frames(child_list(ii,1));
        frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child_list(ii,1),:); ...
        movieInfo.orgCoord(child_list(ii,2),:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
        for jj = 1:length(parent_loc)
            parent_cost(jj) = ovDistanceRegion(movieInfo.vox{detect_divPair(parent_loc(jj),1)},...
            [movieInfo.vox{child_list(ii,1)}; movieInfo.vox{child_list(ii,2)}], frame_shift);
        end
        [~, parent_min] = min(parent_cost);
        parent_remove_flag(parent_loc) = 1;
        parent_remove_flag(parent_loc(parent_min)) = 0;    
    end
end
detect_divPair(logical(parent_remove_flag), :) = [];



%% check division detection accuracy
[TP, TN, FP, gt_divPair, detect_code_list] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) ;
fprintf('TP: %d  TN: %d  FP: %d\n', TP, TN, FP);
% 5: wrong association
% 11: FP detection
% 6 17 32 36 44: postpone 
% 26: sudden appear cell
%% convert FP to movieInfo format for visiualization
tp_flag = ismember(detect_divPair, gt_divPair, 'rows');
fp_pair = detect_divPair(~tp_flag,:);
fp_tracks = cell(size(fp_pair,1), 1);
movieInfo_temp = movieInfo;
movieInfo_temp.parents = cell(size(movieInfo.parents));
movieInfo_temp.kids = cell(size(movieInfo.kids));
for ii = 1:size(fp_pair,1)
    fp_tracks{ii} = fp_pair(ii,:);
    movieInfo_temp.parents{fp_pair(ii,2)} = fp_pair(ii,1);
    movieInfo_temp.kids{fp_pair(ii,1)} = fp_pair(ii,2);
end
movieInfo_temp.tracks = fp_tracks;
mat2tgmm(movieInfo_temp, '/work/Nova/embryo_res_folder/mengfan_data_res/temp');

%%
function z_score = division_gap_test(movieInfo, embryo_vid, child1, child2)
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
        z_score = 100;
    else
        z_score = 0;
    end
end

function [TP, TN, FP, gt_divPair, detect_code_list] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) 
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
    % get postponed gt pairs
    gt_postPair = nan(size(gt_list,1), 4);
    for ii = 1:length(gt_list)
        gt_pair = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_list{ii}),2);
        for jj = 1:3
            if jj == 1
                cell_ind1 = str2num(gt_pair{1});
                cell_ind1 = find(spot_table(:,1) == cell_ind1);
                cell_loc1 = spot_table(cell_ind1, 3:5);
                cell_ind2 = str2num(gt_pair{2});
                cell_ind2 = find(spot_table(:,1) == cell_ind2);
                cell_loc2 = spot_table(cell_ind2, 3:5);
                cell_loc = round((cell_loc1 + cell_loc2)/2);
                gt_postPair(ii, 4) = gt_divPair(ii, 4)+1;
            else
                cell_ind = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_pair{jj-1}),2);
                if length(cell_ind) == 1
                    cell_ind = find(spot_table(:,1) == str2num(cell_ind{1}));
                    cell_loc = round(spot_table(cell_ind, 3:5));
                elseif isempty(cell_ind)
                    cell_loc = nan(1,3);
                else
                    error('Consective division.');
                end
            end
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            gt_postPair(ii,jj) = cell_idx;
        end
    end
    % add missed annotation
    gt_addPair = [53044942 51210387 55802059 40;
                  116506827 113737225 117436109 60];
    gt_divPair = [gt_divPair; gt_addPair];
    gt_addPair(:,4) = gt_addPair(:,4) + 1;
    gt_postPair = [gt_postPair; gt_addPair];
    % remove one incorrect annotation
    error_loc = find(cellfun(@(x)strcmp(x, '34390'), gt_list));
    gt_list(error_loc) = [];
    gt_divPair(error_loc,:) = [];
    gt_postPair(error_loc,:) = [];
    gt_correctPair = [41064118 41989558 41059316 17;];
    gt_divPair = [gt_divPair; gt_correctPair];
    gt_correctPair(:,4) = gt_correctPair(:,4) + 1;
    gt_postPair = [gt_postPair; gt_correctPair];
  

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
        % last try: postpone gt to 1 frame later
        parent_id = findMatchCell(movieInfo, gt_postPair(ii,1), gt_postPair(ii,4));
        detect_loc = find(parent_id == detect_divPair(:,1));
        child1_id = findMatchCell(movieInfo, gt_postPair(ii,2), gt_postPair(ii,4)+1);
        child2_id = findMatchCell(movieInfo, gt_postPair(ii,3), gt_postPair(ii,4)+1);
        if ismember(child1_id, detect_divPair(detect_loc,2:3)) && ...
           ismember(child2_id, detect_divPair(detect_loc,2:3))              % case 1: correct
            gt_code_list(ii) = 1;
            detect_code_list(detect_loc) = 1;
            continue;
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
end
