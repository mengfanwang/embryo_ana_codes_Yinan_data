clc;clear;close all;
% merge the groun truth result to new result

%% get ground truth
gt_path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(gt_path, 0);
% link table: s t
link_table = link_table(:,1:2);
% spot table: id t y x z
% cell fate [NC adaxial] filter lineage [adaxial notochord mixed]
% spot_table = spot_table(:,1:5);

% convert to tracking result corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1/2 1/2 1/5.86] + 1;

%% convert to tracks
d_in = indegree(G);
d_out = outdegree(G);

track_num = sum(d_in == 0) + sum(d_out == 2);
component_id_list = conncomp(G, 'Type', 'weak')';
tracks = cell(track_num,1);
link_flag = zeros(size(link_table,1),1);

track_id = 0;
for component_id = 1:max(component_id_list)
    component_id
    s_ind = find(component_id_list == component_id & d_in == 0);
    s_node = G.Nodes.Name{s_ind};
    t_ind_list = find(component_id_list == component_id & d_out == 0);
    t_node_list = G.Nodes.Name(t_ind_list);
    for ii = 1:length(t_node_list)
        track_id = track_id + 1;
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
                if isempty(tracks{track_id})
                    tracks{track_id} = [parent child];
                else
                    tracks{track_id}(end+1) = child;
                end
                link_flag(link_ind) = 1;
            end
        end
    end
end

% sorting along time 
[~, time_order] = sort(spot_table(:,2));
spot_table = spot_table(time_order,:);

node_num = size(spot_table,1);
node_ind_list = zeros(max(spot_table(:,1)),1);  % node index in G.Nodes
for ii = 1:node_num
   node_ind_list(spot_table(ii,1)) = ii;
end
for track_id = 1:length(tracks)
    for ii = 1:length(tracks{track_id})
        tracks{track_id}(ii) = node_ind_list(tracks{track_id}(ii));
    end
end


%% load movieInfo
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');

%% find missing cell
data_size = [960 960 181];
missing_list = zeros(node_num, 1);
mapping_list = nan(node_num, 1);
fprintf('Check missing cell...\n');
for ii = 1:node_num
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, node_num);
    end
    ii_loc = round(spot_table(ii, 3:5));
    ii_idx = sub2ind(data_size, ii_loc(2), ii_loc(1), ii_loc(3));
    ii_time = spot_table(ii, 2);
    ind_candi = find(movieInfo.frames == ii_time);
    ind_match = find(cellfun(@(x) ismember(ii_idx, x), movieInfo.voxIdx(ind_candi)));
    if isempty(ind_match)
        missing_list(ii) = 1;
    else
        ind_match = ind_match + sum(movieInfo.n_perframe(1:ii_time-1));
        if ismember(ind_match, mapping_list)
            missing_list(ii) = 1;
        else
            mapping_list(ii) = ind_match;
        end
    end
end

movieInfo_temp = movieInfo;
movieInfo = struct();
movieInfo.xCoord = movieInfo_temp.xCoord;
movieInfo.yCoord = movieInfo_temp.yCoord;
movieInfo.zCoord = movieInfo_temp.zCoord;
movieInfo.frames = movieInfo_temp.frames;
movieInfo.n_perframe = movieInfo_temp.n_perframe;
movieInfo.parents = movieInfo_temp.parents;
movieInfo.tracks = movieInfo_temp.tracks;

%% merge gt tracks to movieInfo
% clear; load('/work/Nova/embryo_res_folder/mengfan_data_res/gt2result/temp.mat');
missing_list = find(missing_list);
%movieInfo require: x/y/z coord frames n_perframe parents tracks

t = max(movieInfo.frames);
% add missing cell
for ii = 1:length(missing_list)
    missing_id = missing_list(ii);
    missing_frame = spot_table(missing_id,2);
    % insert to the end of current frame
    id_inMovieInfo = sum(movieInfo.n_perframe(1:missing_frame)) + 1;
    % update movieInfo
    movieInfo.xCoord(id_inMovieInfo+1:end+1) = movieInfo.xCoord(id_inMovieInfo:end);
    movieInfo.xCoord(id_inMovieInfo) = spot_table(missing_id, 3);
    movieInfo.yCoord(id_inMovieInfo+1:end+1) = movieInfo.yCoord(id_inMovieInfo:end);
    movieInfo.yCoord(id_inMovieInfo) = spot_table(missing_id, 4);
    movieInfo.zCoord(id_inMovieInfo+1:end+1) = movieInfo.zCoord(id_inMovieInfo:end);
    movieInfo.zCoord(id_inMovieInfo) = spot_table(missing_id, 5);
    if movieInfo.frames(id_inMovieInfo-1) ~= missing_frame || (missing_frame < t ...
        && movieInfo.frames(id_inMovieInfo) ~= missing_frame + 1)
        error('Error happens!');
    end
    movieInfo.frames(id_inMovieInfo+1:end+1) = movieInfo.frames(id_inMovieInfo:end);
    movieInfo.frames(id_inMovieInfo) = missing_frame;
    movieInfo.n_perframe(missing_frame) = movieInfo.n_perframe(missing_frame) + 1;
    % id after id_inMovieInfo should + 1. Notice that missing cell is not
    % in track
    movieInfo.parents(id_inMovieInfo+1:end+1) = movieInfo.parents(id_inMovieInfo:end);
    movieInfo.parents{id_inMovieInfo} = [];
    for jj = 1:length(movieInfo.parents)
        if movieInfo.parents{jj} >= id_inMovieInfo
            movieInfo.parents{jj} = movieInfo.parents{jj} + 1;
        end
    end
    for jj = 1:length(movieInfo.tracks)
        for kk = 1:length(movieInfo.tracks{jj})
            if movieInfo.tracks{jj}(kk) >= id_inMovieInfo
                movieInfo.tracks{jj}(kk) = movieInfo.tracks{jj}(kk) + 1;
            end
        end
    end
    % update mapping_list
    mapping_list(missing_id) = id_inMovieInfo;
    mapping_list(mapping_list > id_inMovieInfo) = mapping_list(mapping_list > id_inMovieInfo) + 1;
end

% transfer gt tracks to movieInfo index
for ii = 1:length(tracks)   
    tracks{ii} = mapping_list(tracks{ii});
end

% update anootated cells locations
for ii = 1:size(spot_table,1)
    movieInfo.xCoord(mapping_list(ii)) = spot_table(ii,3);
    movieInfo.yCoord(mapping_list(ii)) = spot_table(ii,4);
    movieInfo.zCoord(mapping_list(ii)) = spot_table(ii,5);
end

%% merge tracks
tracks_gt = tracks;
tracks_cur = movieInfo.tracks;
tracks_candi = {};
for ii = 1:length(tracks_cur)
    collision_flag = false(size(tracks_cur{ii}));
    for jj = 1:length(tracks_gt)
        collision_flag = collision_flag | ismember(tracks_cur{ii}, tracks_gt{jj});
    end
    if any(collision_flag)
        survive_loc = find(~collision_flag);
        if ~isempty(survive_loc)
            break_loc = find(diff(survive_loc)>1);
            if isempty(break_loc)
                tracks_candi{end+1} = tracks_cur{ii}(survive_loc);
            else
                tracks_candi{end+1} = tracks_cur{ii}(survive_loc(1):survive_loc(break_loc(1)));
                for jj = 1:length(break_loc)-1
                    tracks_candi{end+1} = tracks_cur{ii}(survive_loc(break_loc(jj)+1):survive_loc(break_loc(jj+1)));
                end
                tracks_candi{end+1} = tracks_cur{ii}(survive_loc(break_loc(end)+1):survive_loc(end));
            end
        end
        tracks_cur{ii} = [];
    end
end
tracks_cur(cellfun(@isempty,tracks_cur)) = [];
movieInfo.tracks = [tracks_cur; tracks_candi'; tracks_gt;];
% update parent
movieInfo.parents = cell(size(movieInfo.parents));
for ii = 1:length(movieInfo.tracks)
    for jj = 2:length(movieInfo.tracks{ii})
        movieInfo.parents{movieInfo.tracks{ii}(jj)} = movieInfo.tracks{ii}(jj-1);
    end
end

%% add tags in the form of radius
radius = [12.5 12 11.5 11 10.5];  % cell fate [1 2]  lineage [1 2 3]
% cell fate has higher priority
category = zeros(size(movieInfo.xCoord));
for ii = 1:node_num
    if spot_table(ii,6) == 1
        category(mapping_list(ii)) = 1;
    elseif spot_table(ii,7) == 1
        category(mapping_list(ii)) = 2;
    elseif spot_table(ii,9) == 1
        category(mapping_list(ii)) = 3;
    elseif spot_table(ii,10) == 1
        category(mapping_list(ii)) = 4;
    elseif spot_table(ii,11) == 1
        category(mapping_list(ii)) = 5;
    end
end

mat2tgmm_tags(movieInfo, '/work/Nova/embryo_res_folder/mengfan_data_res/gt2result/tgmm', radius, category);

%% check all annotated cells are in tracks
inTrack_flag = zeros(node_num,1);
for ii = 1:length(movieInfo.tracks)
    inTrack_flag(ismember(mapping_list, movieInfo.tracks{ii})) = 1;
end

%% check no linked cells
single_end = [];
for ii = 1:node_num
    if ~ismember(spot_table(ii,1), link_table(:,1)) && ~ismember(spot_table(ii,1), link_table(:,2))
        if isempty(movieInfo.parents{mapping_list(ii)})
            single_end = [single_end ii];
        end
    end
end

%% check frames
borken_track = [];
frames = cell(size(movieInfo.tracks));
for ii = 1:length(frames)
    frames{ii} = movieInfo.frames(movieInfo.tracks{ii});
    if any(diff(diff(frames{ii}))~=0)
        borken_track = [borken_track ii];
    end
end

%% plot ground truth only
% movieInfo require: x/y/z coord n_perframe parents tracks
% movieInfo_gt = struct();
% movieInfo_gt.tracks = tracks;
% movieInfo_gt.xCoord = spot_table(:,3);
% movieInfo_gt.yCoord = spot_table(:,4);
% movieInfo_gt.zCoord = spot_table(:,5);
% movieInfo_gt.n_perframe = zeros(max(spot_table(:,2)),1);
% for tt = 1:max(spot_table(:,2))
%     movieInfo_gt.n_perframe(tt) = sum(spot_table(:,2) == tt);
% end
% movieInfo_gt.parents = cell(node_num,1); 
% for track_id = 1:length(tracks)
%     for ii = 2:length(tracks{track_id})
%         movieInfo_gt.parents{tracks{track_id}(ii)} = tracks{track_id}(ii-1);
%     end
% end

% mat2tgmm_tags(movieInfo_gt, '/work/Nova/embryo_res_folder/mengfan_data_res/gt_temp/', radius, table);
