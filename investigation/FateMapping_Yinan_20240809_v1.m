clc;clear;close all;
% stitch tracking segments 
addpath ../src_code_matlab/ ../TGMM_wrapper/

%% parameter setting
files = {{'000', '050'}, {'041', '090'}, {'081', '133'}};

max_dist = 30;
resolution = [2 2 6];
size_ratio = 1;
max_nei = 2;
tp_margin = 5; % overlap/2. don't use the last several files 

tps = str2num(files{1}{1}):str2num(files{end}{2});
pre_tps = -1;
    
path = '/work/Mengfan/Embryo/24-08-09_Yinan_v1';
data_folder = fullfile(path, 'data');
file_path = fullfile(path, 'Tracking');
translation_path = fullfile(path, 'Registration');
[driftInfo, drift] = loadRegInfo(translation_path, tps(1:end-1));

%%
% ss = 1;
% m1 = load(fullfile(file_path, [files{ss}{1} '_' files{ss}{2}], 'movieInfo.mat'));
% m1 = m1.movieInfo;
% m2 = load(fullfile(file_path, [files{ss+1}{1} '_' files{ss+1}{2}], 'movieInfo.mat'));
% m2 = m2.movieInfo;

%% stitch part 1: merge all files into one

movieInfoAll.n_perframe = [];
movieInfoAll.frames = [];
movieInfoAll.vox = {};
movieInfoAll.voxIdx = {};
movieInfoAll.tracks = {};
movieInfoAll.orgCoord = [];
movieInfoAll.xCoord = [];
movieInfoAll.yCoord = [];
movieInfoAll.zCoord = [];
movieInfoAll.parents = {};  
movieInfoAll.kids = {};

for ff = 1:length(files)
    ff
    if ff == length(files)
        tp_margin = 0;
    end
load(fullfile(file_path, [files{ff}{1} '_' files{ff}{2}], 'movieInfo.mat'));

% check tps to process in cur file
tps_file = str2num(files{ff}{1}):str2num(files{ff}{2});
if ~ismember(pre_tps+1, tps_file)
    error('Files are not overlapped!');
end
frame_range = find(tps_file(1:end-tp_margin) > pre_tps);
n_cumsum = [0; cumsum(movieInfo.n_perframe)];

% adjust parent kids tracks relationship
track_remove_flag = zeros(length(movieInfo.tracks),1);
for ii = 1:length(movieInfo.tracks)
    remove_flag = zeros(length(movieInfo.tracks{ii}),1);
    for jj = 1:length(movieInfo.tracks{ii})-1
        parent = movieInfo.tracks{ii}(jj);
        kid = movieInfo.tracks{ii}(jj+1);
        if movieInfo.frames(parent) < frame_range(1)
            remove_flag(jj) = 1;
            movieInfo.parents{kid} = [];
            movieInfo.kids{parent} = [];
        end
    end
    for jj = length(movieInfo.tracks{ii}):-1:2
        parent = movieInfo.tracks{ii}(jj-1);
        kid = movieInfo.tracks{ii}(jj);
        if movieInfo.frames(kid) > frame_range(end)
            remove_flag(jj) = 1;
            movieInfo.parents{kid} = [];
            movieInfo.kids{parent} = [];
        end
    end
    movieInfo.tracks{ii}(logical(remove_flag)) = [];
    if isempty(movieInfo.tracks{ii})
        track_remove_flag(ii) = 1;
    end
end
movieInfo.tracks(logical(track_remove_flag)) = [];


% key info: n_perframe frames vox voxIdx tracks orgCoord x/y/zCoord parents kids
cell_range = n_cumsum(frame_range(1))+1:n_cumsum(frame_range(end)+1);
% correct cell index
ind_correction = -n_cumsum(frame_range(1)) + length(movieInfoAll.voxIdx);
for ii = 1:length(movieInfo.tracks)
    movieInfo.tracks{ii} = movieInfo.tracks{ii} + ind_correction;
end
for ii = 1:length(movieInfo.voxIdx)
    if ~isempty(movieInfo.parents{ii})
        movieInfo.parents{ii} = movieInfo.parents{ii} + ind_correction;
    end
    if ~isempty(movieInfo.kids{ii})
        movieInfo.kids{ii} = movieInfo.kids{ii} + ind_correction;
    end
end
movieInfoAll.frames = [movieInfoAll.frames; movieInfo.frames(cell_range(1):cell_range(end)) - ...
                                            frame_range(1) + length(movieInfoAll.n_perframe) + 1];
movieInfoAll.n_perframe = [movieInfoAll.n_perframe; movieInfo.n_perframe(frame_range(1):frame_range(end))];
movieInfoAll.vox = [movieInfoAll.vox; movieInfo.vox(cell_range(1):cell_range(end))];
movieInfoAll.voxIdx = [movieInfoAll.voxIdx; movieInfo.voxIdx(cell_range(1):cell_range(end))];
movieInfoAll.tracks = [movieInfoAll.tracks; movieInfo.tracks];
movieInfoAll.orgCoord = [movieInfoAll.orgCoord; movieInfo.orgCoord(cell_range(1):cell_range(end), :)];
movieInfoAll.xCoord = [movieInfoAll.xCoord; movieInfo.xCoord(cell_range(1):cell_range(end))];
movieInfoAll.yCoord = [movieInfoAll.yCoord; movieInfo.yCoord(cell_range(1):cell_range(end))];
movieInfoAll.zCoord = [movieInfoAll.zCoord; movieInfo.zCoord(cell_range(1):cell_range(end))];
movieInfoAll.parents = [movieInfoAll.parents; movieInfo.parents(cell_range(1):cell_range(end))];
movieInfoAll.kids = [movieInfoAll.kids; movieInfo.kids(cell_range(1):cell_range(end))];
pre_tps = pre_tps + length(frame_range);
end
%% stitch part 2: try the best to find parents for every cell
clear movieInfo
n_cumsum = [0; cumsum(movieInfoAll.n_perframe)];
assigned = zeros(length(movieInfoAll.n_perframe),1);
for ii = 1:length(movieInfoAll.tracks)
    if mod(ii, 1000) == 0
        fprintf("%d/%d\n", ii, length(movieInfoAll.tracks));
    end
    if isempty(movieInfoAll.tracks{ii})
        continue
    end
    track_head = movieInfoAll.tracks{ii}(1);
    if movieInfoAll.frames(track_head) == 1
        continue
    end
    tt = movieInfoAll.frames(track_head);
    cell_loc = movieInfoAll.orgCoord(track_head,:);
    cell_loc = mapCellLoc(cell_loc, tt-1, tt, driftInfo, drift);
    dist_pair = pdist2(movieInfoAll.orgCoord(n_cumsum(tt-1)+1:n_cumsum(tt), :).*resolution, cell_loc.*resolution);
    
    % get neighbors and neighbor distance
    [neighbor_dist, neighbor_candidate] = sort(dist_pair);
    neighbor_dist = neighbor_dist(1:min(max_nei, length(dist_pair)));
    neighbor_candidate = neighbor_candidate(1:min(max_nei, length(dist_pair)));
    neighbors = neighbor_candidate(neighbor_dist < max_dist);
    neighbors = neighbors + n_cumsum(tt-1);
    dist2nei = neighbor_dist(neighbor_dist < max_dist);

    % the neighbirs within distance and k-nearest are candidates
    % then check them according to order]
    for jj = 1:length(neighbors)
        parent = neighbors(jj);
        parent_kids = movieInfoAll.kids{parent};
        kid_num = length(parent_kids);
        % if kids < 2 and one of the size < 0.8 -> connect
        if kid_num == 0
            movieInfoAll.tracks{ii} = [parent; movieInfoAll.tracks{ii}];
            movieInfoAll.parents{track_head} = parent;
            movieInfoAll.kids{parent} = [movieInfoAll.kids{parent}; track_head];
            assigned(tt) = assigned(tt)+1;
            break
        end
    end
end


%% writing
mkdir(fullfile(file_path, 'stitch'));
save(fullfile(file_path, 'stitch', 'movieInfo.mat'), 'movieInfoAll','-v7.3');
mastodon_dir = fullfile(file_path, 'stitch', 'mastodon');
mkdir(mastodon_dir);
mat2tgmm(movieInfoAll, fullfile(mastodon_dir, 'tgmm_format'));
tif2bdv(data_folder, fullfile(mastodon_dir, 'embryo_data_h5'), generate_tps_str(tps), [], []);


function cell_loc_new = mapCellLoc(cell_loc, start_tp, end_tp, driftInfo, drift)
    cell_loc_new = cell_loc;
    for tt = end_tp:-1:start_tp
        y_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.y_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        x_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.x_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        z_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.z_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        cell_loc_new = cell_loc_new + [y_bias x_bias z_bias];
    end
end

function [driftInfo, drift] = loadRegInfo(translation_path, tps)
    
    timepts_to_process = generate_tps_str(tps);
    driftInfo.grid_size = 32;    % vector numbers each dim, currently cube only
    driftInfo.batch_size = [30 30 8];
    
    [y_batch, x_batch, z_batch] = meshgrid(0:driftInfo.grid_size+1);
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2;
    
    driftInfo.y_batch = y_batch;
    driftInfo.x_batch = x_batch;
    driftInfo.z_batch = z_batch;  % vector locations in original image with one padding
    
    grid_size = driftInfo.grid_size;
    drift.x_grid = cell(length(timepts_to_process),1);
    drift.y_grid = cell(length(timepts_to_process),1);
    drift.z_grid = cell(length(timepts_to_process),1);
    for ii = 1:length(timepts_to_process)
        load(fullfile(translation_path, timepts_to_process(ii)+'.mat'), 'phi_current_vec');
        x_grid = padarray(reshape(phi_current_vec(1:3:end-2), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        y_grid = padarray(reshape(phi_current_vec(2:3:end-1), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        z_grid = padarray(reshape(phi_current_vec(3:3:end), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
        drift.x_grid{ii} = x_grid;
        drift.y_grid{ii} = y_grid;
        drift.z_grid{ii} = z_grid;
    end
end