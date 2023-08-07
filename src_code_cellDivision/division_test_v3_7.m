% clc;clear;close all;
% load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0504_000_191/movieInfo.mat');

% use v3.4 pipeline
% process the whole data to reduce FP
clearvars -except movieInfo embryo_vid
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

%% find all detections and candidate divisions in movieInfo
t = 192;
data_size = [960 960 181];
im_resolution = [1 1 5.86/2];
Ci = abs(movieInfo.Ci(1));
cell_num = length(movieInfo.xCoord);
div_candiList = cell(cell_num,1);
div_candiList_2nd = cell(cell_num,1);
fprintf('Find cells and candidate divisions in movieInfo...\n'); 
% get heads of tracks
head_list = nan(length(movieInfo.tracks),1);
head_list_2nd = nan(length(movieInfo.tracks),1);
length_thre = 10;
for ii = 1:length(movieInfo.tracks)
    %
    track = movieInfo.tracks{ii};
    if length(track) >= length_thre
        head_list(ii) = track(1);
        z_inten = zeros(length_thre,1);
        mean_inten = zeros(length_thre,1);
        for jj = 1:length_thre
            voxInten = embryo_vid{movieInfo.frames(track(jj))}(movieInfo.voxIdx{track(jj)});
            mean_inten(jj) = mean(voxInten);
            [vidOut, vIdx, ~, ~] = crop3D(embryo_vid{movieInfo.frames(track(jj))},...
                movieInfo.voxIdx{track(jj)}, [4 4 2]);
            vBin = zeros(size(vidOut));
            vBin(ismember(vIdx, movieInfo.voxIdx{track(jj)})) = 1;
        
            z_inten(jj) = bgZtest(vidOut, vBin, 3);
        end
        if z_inten(1) == min(z_inten) && mean_inten(1) == min(mean_inten)
            head_list_2nd(ii) = track(2);
        end
    end
end
tic;
for parent_id = 1:cell_num
    if mod(parent_id,1000) == 0
        fprintf('%d / %d\n', parent_id, cell_num);
    end
    parent_time = movieInfo.frames(parent_id);
    % if have child already, it should be small
    child_small_flag = true;
    child_id = movieInfo.kids{parent_id};
    if ~isempty(child_id)
        if length(movieInfo.voxIdx{child_id}) >= ...
            length(movieInfo.voxIdx{parent_id}) * 0.8   
            child_small_flag = false;
        end
%         if sum(embryo_vid{parent_time+1}(movieInfo.voxIdx{child_id})) >= ...
%             sum(embryo_vid{parent_time}(movieInfo.voxIdx{parent_id})) * 0.8
%             child_small_flag = true;
%         end
    end
    if parent_time ~= t && child_small_flag
        % find the 5-nearest neighbor
        div_candi = findNeighbor(movieInfo, parent_id, im_resolution);
        % only heads are included   
        div_candi1 = div_candi(ismember(div_candi,head_list));
        div_candiList{parent_id} = div_candi1;

        % second head included
        div_candi2 = div_candi(ismember(div_candi,head_list_2nd));
        div_candiList_2nd{parent_id} = div_candi2;
    end
end
toc                                    

%% convert candiList to division cases
div_num = 0;
detect_divPair = zeros(1000000,3);
div_2ndnum = 0;    % temporary for test
detect_2ndPair = zeros(10000,3);
for ii = 1:cell_num
    if mod(ii,10000) == 0
        fprintf('%d / %d\n', ii, cell_num);
    end
    parent_id = ii;
    detect_flag = false;
    if ~isempty(movieInfo.kids{parent_id}) && ...
            ~isempty(div_candiList{ii})         % case1: one child exists
        child1_id = movieInfo.kids{parent_id};
        children_pair = [repmat(child1_id, [length(div_candiList{ii}) 1])...
                             div_candiList{ii}'];
        child_time = movieInfo.frames(child1_id);
        detect_divPair_tmp = candidateFiltering(movieInfo, embryo_vid{child_time}, parent_id, ...
                                children_pair, 1, data_size, im_resolution);  
        if ~isempty(detect_divPair_tmp)
            div_num = div_num + 1;
            detect_divPair(div_num,:) = detect_divPair_tmp;
            detect_flag = true;
        end
    elseif length(div_candiList{ii}) >= 2       % case2: no child exists
        children_pair = nchoosek(div_candiList{ii}, 2);
        
        detect_divPair_tmp = candidateFiltering(movieInfo, embryo_vid{child_time}, parent_id, ...
                                children_pair, 1, data_size, im_resolution);   
        if ~isempty(detect_divPair_tmp)
            div_num = div_num + 1;
            detect_divPair(div_num,:) = detect_divPair_tmp;
            detect_flag = true;
        end
    end
    if ~detect_flag
        % try to find the second head
        if ~isempty(movieInfo.kids{parent_id})  && ...
                ~isempty(div_candiList_2nd{parent_id})
            % repeat the same step
            child1_id = movieInfo.kids{parent_id};
            children_pair = [repmat(child1_id, [length(div_candiList_2nd{parent_id}) 1])...
                                div_candiList_2nd{parent_id}'];
            child_time = movieInfo.frames(child1_id);
            detect_divPair_tmp = candidateFiltering(movieInfo, embryo_vid{child_time}, parent_id, ...
                                children_pair, 1, data_size, im_resolution); 
            if ~isempty(detect_divPair_tmp)
                div_num = div_num + 1;
                detect_divPair(div_num,:) = detect_divPair_tmp;
            end
        end
    end
end
detect_divPair = detect_divPair(1:div_num,:);
detect_divPair = unique(detect_divPair, 'rows');

detect_2ndPair = detect_2ndPair(1:div_2ndnum, :);

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
        dist2parent(jj) = norm( (mean(movieInfo.vox{parent_id}) + frame_shift -...
            mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}])) .* im_resolution);
    end
     [~, minIdx] = min(dist2parent);
    parent_remove_flag(parent_loc) = 1;
    parent_remove_flag(parent_loc(minIdx)) = 0;    

end
detect_divPair(logical(parent_remove_flag), :) = [];

%% keep annotated cell only for evaluation
cell_match_list = nan(size(spot_table,1),1);
for ii = 1:size(spot_table,1)
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, size(spot_table,1));
    end
    cell_loc = round(spot_table(ii, 3:5));
    cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
    cell_time = spot_table(ii, 2);
    ind_candi = find(movieInfo.frames == cell_time);
    ind_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(ind_candi)));
    if ~isempty(ind_match)
        detection_id = ind_match + sum(movieInfo.n_perframe(1:cell_time-1));
        cell_match_list(ii) = detection_id;
    end
end
%%
in_gt_flag = false(size(detect_divPair,1),1);
for ii = 1:size(detect_divPair,1)
    if ismember(detect_divPair(ii,1), cell_match_list)
        in_gt_flag(ii) = true;
    end
end
detect_divPair = detect_divPair(in_gt_flag,:);
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
[TP, FN, FP, gt_divPair, gt_code_list, detect_code_list] = ...
    div_acc(detect_allPair, G, movieInfo, spot_table, data_size);
fprintf('TP: %d  FN: %d  FP: %d\n', TP, FN, FP);

%% FN investigation
gt_fnPair = gt_divPair(gt_code_list~=1,:);
cell_match_flag = nan(size(gt_fnPair,1), 3);
size_ratio = nan(size(gt_fnPair,1), 2);
track_loc = nan(size(gt_fnPair,1), 2);
track_length = nan(size(gt_fnPair,1), 2);
parent_track_flag = false(size(gt_fnPair,1), 2);
not_neighbor_flag = false(size(gt_fnPair,1), 2);
too_large_flag = false(size(gt_fnPair,1), 1);
not_nearest_parent_flag = false(size(gt_fnPair,1), 2);
no_gap_flag = false(size(gt_fnPair,1), 1);
dist2parent = inf(size(gt_fnPair,1), 2);
dist2nearest = inf(size(gt_fnPair,1), 2);
for ii = 1:size(gt_fnPair,1)
    parent_id = findMatchCell(movieInfo, gt_fnPair(ii,1), gt_fnPair(ii,4));
    child1_id = findMatchCell(movieInfo, gt_fnPair(ii,2), gt_fnPair(ii,4)+1);
    child2_id = findMatchCell(movieInfo, gt_fnPair(ii,3), gt_fnPair(ii,4)+1);
    cell_match_flag(ii,:) = [parent_id child1_id child2_id];
    if ~isnan(parent_id)
        parent_track = movieInfo.particle2track(parent_id,1);
        neighbors = findNeighbor(movieInfo, parent_id, im_resolution);
        if ~isempty(movieInfo.kids{parent_id})
            if length(movieInfo.voxIdx{movieInfo.kids{parent_id}}) > ...
                length(movieInfo.voxIdx{parent_id}) * 0.8
                too_large_flag(ii) = true;
            end
        end
        if ~isnan(child1_id) && ~isnan(child2_id)
            no_gap_flag(ii) = division_gap_test(movieInfo, zeros(data_size), child1_id, child2_id);
            no_gap_flag(ii) = ~no_gap_flag(ii);
        end
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
            end
        end
    end
end

fprintf('Missing cell #:%d %f\n', sum(any(isnan(cell_match_flag)')), ...
    sum(any(isnan(cell_match_flag)'))/size(gt_fnPair,1));

% too_large_flag = size_ratio(:,1)>0.8;
fprintf('Large child #:%d %f\n', sum(too_large_flag), ...
    sum(too_large_flag)/size(gt_fnPair,1));

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

fprintf('Gap test not pass#:%d %f\n', sum(no_gap_flag), ...
    sum(no_gap_flag)/size(gt_fnPair,1));

sum(over_merge_flag | any(too_large_flag, 2) | any(isnan(cell_match_flag),2) ...
    | any(not_head_flag,2) | any(too_short_flag,2) | any(not_neighbor_flag,2) | sum(no_gap_flag)) 
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

function detect_divPair = candidateFiltering(movieInfo, embryo_vid, parent_id, ...
                            children_pair, gap_test_flag, data_size, im_resolution)  
% should statisify the following creteria:
% 1.parent must be the nearest one to the center of the pair of children
% 2.children must split (depends)
% 3.two child size ratio < 3
    child_time = movieInfo.frames(parent_id)+1;
    dist2parent = inf(size(children_pair,1),1);
    for jj = 1:size(children_pair,1)
        child1_id = children_pair(jj,1);
        child2_id = children_pair(jj,2);
        frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
            movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
        dist2parent(jj) = norm( (mean(movieInfo.vox{parent_id}) + frame_shift -...
            mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}])) .* im_resolution);
    end
    [~, minIdx] = min(dist2parent);
    child1_id = children_pair(minIdx, 1);
    child2_id = children_pair(minIdx, 2);

    if gap_test_flag
        gap_flag = division_gap_test(movieInfo, zeros(data_size), child1_id, child2_id);
    else
        gap_flag = true;
    end
    child_size = [length(movieInfo.voxIdx{child1_id})...
                  length(movieInfo.voxIdx{child2_id})];
    size_flag = (max(child_size) / min(child_size)) < 3;
    
    Ci = 33.7382;
    phat = [2.0557    0.7412];
    frame_shift1 = getNonRigidDrift([0 0 0], movieInfo.orgCoord(child1_id,:),...
        child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
    [ov_dist1, ~, ~] = ovDistanceRegion(movieInfo.vox{parent_id},...
        movieInfo.vox{child1_id}, frame_shift1);
    ov_dist1 = overlap2cost(ov_dist1, phat);
    frame_shift2 = getNonRigidDrift([0 0 0], movieInfo.orgCoord(child2_id,:),...
        child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
    [ov_dist2, ~, ~] = ovDistanceRegion(movieInfo.vox{parent_id},...
        movieInfo.vox{child2_id}, frame_shift2);
    ov_dist2 = overlap2cost(ov_dist2, phat);
    dist_flag = (ov_dist1 < Ci) & (ov_dist2 < Ci);

%     im1 = embryo_vid(movieInfo.voxIdx{child1_id});
%     im2 = embryo_vid(movieInfo.voxIdx{child2_id});
%     [~, p] = ttest2(im1,im2);
%     child_diff_flag = p > normcdf(-5)*2;  %z-score = 5

    if gap_flag && size_flag %&& child_diff_flag
        detect_divPair = [parent_id children_pair(minIdx,:)];
    else
        detect_divPair = [];
    end

end

function gap_flag = division_gap_test(movieInfo, embryo_vid, child1, child2)
    % get local region
    [vidOut, vIdx, ~, ~] = crop3D(embryo_vid, [movieInfo.voxIdx{child1}; movieInfo.voxIdx{child2}], [0 0 0]);
    vBin = zeros(size(vidOut));
    vBin(ismember(vIdx, movieInfo.voxIdx{child1})) = 1;
    vBin(ismember(vIdx, movieInfo.voxIdx{child2})) = 1;

    cc = bwconncomp(vBin, 6);
    if cc.NumObjects == 2
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
    error_loc = [find(cellfun(@(x)strcmp(x, '34390'), gt_list));
                 find(cellfun(@(x)strcmp(x, '186025'), gt_list));];
    gt_list(error_loc) = [];
    gt_divPair(error_loc,:) = [];
    gt_correctPair = [41064118 41989558 41059316 17;
                      123787601 123791443 122870796 60];
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


function neighbors = findNeighbor(movieInfo, cur_i, im_resolution)
    % find the nearest neighbors in the next frame with given conditions

    max_dist = 50/2;
    max_nei = 5;
    
    tt = movieInfo.frames(cur_i);
    curCentroid = movieInfo.orgCoord(cur_i,:);
    nextCentroid = movieInfo.orgCoord(movieInfo.frames==(tt+1),:);
    
    % candidate neighbors roughly selection
    drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt, tt+1, movieInfo.drift, movieInfo.driftInfo);
    dist_pair = pdist2((curCentroid + drift).*im_resolution, nextCentroid.*im_resolution); % not accurate but enough
    [neighbor_dist, neighbor_candidate] = sort(dist_pair);
    neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
    neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
    neighbors = neighbor_candidate(neighbor_dist < max_dist);
    neighbors = neighbors + sum(movieInfo.n_perframe(1:tt));
end

function z_score = bgZtest(data, fMap, strel_rad)
    strel_ker = ones(strel_rad*2+1, strel_rad*2+1,ceil(strel_rad/3)*2+1);
    [xx,yy,zz] = ind2sub_direct(size(strel_ker), find(strel_ker));
    dist = sqrt( (xx - strel_rad-1).^2 + (yy - strel_rad-1).^2 + ( (zz - ceil(strel_rad/3)-1)*5 ).^2 );
    strel_ker(dist>=strel_rad) = 0;
    strel_ker = strel(strel_ker);

    fg = data(logical(fMap));
    fMap = imdilate(fMap,strel_ker) - fMap;
    bg = data(logical(fMap));
%     z_score = (mean(fg) - mean(bg))/std(bg)*sqrt(length(fg));
    [mu, sigma] = ordStatApproxKsecWith0s_mat(fg, bg, []);
    L = mean(fg) - mean(bg);
    z_score = (L*sqrt(length(fg)+length(bg)) - mu) / sigma;
end

% function neighbors = findParent(movieInfo, cur_i)
%     % find the nearest neighbors in the next frame with given conditions
% 
% %     max_dist = 50/2;
% %     max_nei = 5;
% %     resolution = [1 1 5.86/2];
%     
%     tt = movieInfo.frames(cur_i);
%     curCentroid = movieInfo.orgCoord(cur_i,:);
%     prevCentroid = movieInfo.orgCoord(movieInfo.frames==(tt-1),:);
%     
%     % candidate neighbors roughly selection
%     drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt-1, tt, movieInfo.drift, movieInfo.driftInfo);
%     dist_pair = pdist2((curCentroid + drift).*resolution, prevCentroid.*resolution); % not accurate but enough
%     [neighbor_dist, neighbor_candidate] = sort(dist_pair);
%     neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
%     neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
%     neighbors = neighbor_candidate(neighbor_dist < max_dist);
%     neighbors = neighbors + sum(movieInfo.n_perframe(1:tt-2));
% end