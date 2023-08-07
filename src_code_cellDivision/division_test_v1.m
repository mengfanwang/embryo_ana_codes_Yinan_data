clc;clear;close all;
% criteria:  relative size < 0.8
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
    if ~isempty(movieInfo.tracks{ii})
        head_list(ii) = movieInfo.tracks{ii}(1);
    end
end
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
        if cell_time ~= t
            div_candiList{ii} = findNeighbor(movieInfo, detection_id)...
                + sum(movieInfo.n_perframe(1:cell_time));
            % only heads are included   
            head_flag = ismember(div_candiList{ii},head_list);
            div_candiList{ii} = div_candiList{ii}(head_flag);
        end
    end
end



%% check division detection accuracy

% div_candiList -> detect_divPair
detect_divPair = zeros(cell_num*10,2);
div_ind = 0;
for ii = 1:cell_num
    for jj = 1:length(div_candiList{ii})
        if length(movieInfo.voxIdx{div_candiList{ii}(jj)}) < ...
                length(movieInfo.voxIdx{detection_list(ii)}) * 0.8
            div_ind = div_ind + 1;
            detect_divPair(div_ind,:) = [detection_list(ii), div_candiList{ii}(jj)];
        end
    end
end
detect_divPair = detect_divPair(1:div_ind,:);
[TP, TN, FP, gt_divPair, detect_divPair] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) ;
fprintf('TP: %d  TN: %d  FP: %d\n', TP, TN, FP);

%% plot the intensity variation
detect_inten = nan(size(detect_divPair,1), 400);   % parent at 100
tp_flag = ismember(detect_divPair, gt_divPair, 'rows');
% get track and it's intensity
cell_inten = cat(1, null_intensity{:});
for ii = 1:size(detect_divPair,1)
    parent = detect_divPair(ii,1);
    child = detect_divPair(ii,2);
    track_id = find(child == head_list);
    if isempty(track_id)
        % already detect before
        for jj = 1:length(movieInfo.tracks)
            if all(ismember(detect_divPair(ii,:), movieInfo.tracks{jj}))
                track = movieInfo.tracks{jj};
                parent_loc = find(track == parent);
                detect_inten_temp = cell_inten(track) ./ null_intensity_mean(movieInfo.frames(track));
                detect_inten(ii, 201-parent_loc: 200-parent_loc+length(track)) ...
                    = detect_inten_temp / mean(detect_inten_temp);
                break;
            end
        end
    else
        track = movieInfo.tracks{track_id};
        detect_inten_temp = cell_inten(track) ./ null_intensity_mean(movieInfo.frames(track));
        detect_inten(ii, 200: 199+length(track)) ...
            = detect_inten_temp / mean(detect_inten_temp);
    end
end
%%
figure(1);hold on
for ii = 1:size(detect_divPair,1)
    ii
    if tp_flag(ii)
        plot(-199:200, detect_inten(ii,:), 'b.');
    else
        plot(-199:200, detect_inten(ii,:), 'r.');
    end
end
%%
figure(2);
fp_child = unique(detect_divPair(~tp_flag,2));
fp_length_flag = zeros(size(fp_child));
for ii = 1:length(fp_child)
    track_id = find(fp_child(ii) == head_list);
    if ~isempty(track_id) && length(movieInfo.tracks{track_id}) >= 20
        fp_length_flag(ii) = 1;
        track = movieInfo.tracks{track_id};
        detect_inten_temp = cell_inten(track) ./ null_intensity_mean(movieInfo.frames(track));
        detect_inten_temp = detect_inten_temp / mean(detect_inten_temp);
        plot(detect_inten_temp, 'r');
        fprintf('ii: %d time: %d', ii, movieInfo.frames(track(1))-1);
        movieInfo.orgCoord(track(1),:).*[2 2 1]
        a = 1;
    end
end
% fp_child(~fp_length_flag) = [];

%%
function [TP, TN, FP, gt_divPair, detect_divPair] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) 
    
    d_out = outdegree(G);
    gt_list = G.Nodes.Name(find(d_out == 2));
    gt_pair = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_list),:);
    gt_divPair = nan(size(gt_pair));
    for ii = 1:size(gt_pair,1)
        for jj = 1:2
            cell_ind = str2num(gt_pair{ii,jj});
            cell_ind = find(spot_table(:,1) == cell_ind);
            cell_loc = round(spot_table(cell_ind, 3:5));
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            cell_time = spot_table(cell_ind, 2);
            cell_candi = find(movieInfo.frames == cell_time);
            cell_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(cell_candi)));
            if isempty(cell_match)
                cell_match = nan;
            else
                cell_match = cell_match + sum(movieInfo.n_perframe(1:cell_time-1));
            end
            if jj == 1
                parent = cell_match;
            else
                child = cell_match;
            end
        end
        gt_divPair(ii,:) = [parent child];
        % add already detected divisions
        if ~isnan(parent)
            if ~isempty(movieInfo.kids{parent})
                detect_divPair(end+1,:) = [parent movieInfo.kids{parent}];
            end
        end
    end
    detect_divPair = unique(detect_divPair,'rows');
    
    detect_flag = ismember(gt_divPair, detect_divPair, 'rows');
    TP = sum(detect_flag);
    TN = length(detect_flag) - TP;
    FP = size(detect_divPair,1) - TP;
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

