clc;clear;close all;
% Try to map every cell into the next time point. If the center is mapped
% into a cell, they are matched.
% Otherwise, we will try to find the appropriate candidates in the
% neighborhoods. (at least one size_ratio < 0.8)
% If there are more than two candidates,  we can treat it as a division.
dbstop if error 

%% load registration info
start_tp = 1;
end_tp = 400;
resolution = [1 1 8];
max_dist = 50;
size_ratio_thre = 0.8;


%% load initial position
load('/work/public/Embryo/Amat2014/Detection/Wei_refine_res/00000.mat');
[cell_loc, cell_size] = refineRes2cellLoc(refine_res);
boundary_label = cell_loc(:,3) >= 115;
cell_loc_all = cell(400,1);
cell_loc_all{1} = [cell_loc boundary_label];

n_perframe = zeros(end_tp,1);
n_perframe(1) = size(cell_loc,1);
n_cumsum = zeros(end_tp+1,1);
n_cumsum(2) = size(cell_loc,1);

frames = zeros(8e6,1);
frames(1:n_cumsum(2)) = 1;
orgCoord = zeros(8e6,3);
orgCoord(1:n_cumsum(2),:) = cell_loc .* [2 2 1];


% %% plot initial positions
% figure(1);title('T = 0');hold on;
% scatter3(cell_loc(~boundary_label,1), cell_loc(~boundary_label,2), cell_loc(~boundary_label,3).*8, 36, [0.5 0.5 0.5], '.');
% scatter3(cell_loc(boundary_label,1), cell_loc(boundary_label,2), cell_loc(boundary_label,3).*8, 72, [0 0.4470 0.7410], '.');
% cell_loc_0 = [cell_loc boundary_label];



timepts_to_process = generate_tps_str(start_tp:end_tp-1);
st_loc = [];
sz_crop = [];
translation_path = '/work/public/Embryo/Amat2014/Registration_reverse';
driftInfo.grid_size = 32;    % vector numbers each dim, currently cube only
driftInfo.batch_size = [30 30 8];

[y_batch, x_batch, z_batch] = meshgrid(0:driftInfo.grid_size+1);
if ~isempty(st_loc)  % adjust if crop
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2 - st_loc(2)/sc_f;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2 - st_loc(1)/sc_f;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2 - st_loc(3);
else
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2;
end
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

%% mapping
% there are two cell locs:
% cell_loc: current cell location from detection
% cell_loc_map : previous cell location mapped to cur
% also two labels:
% boundary_label: label of the prev time point
% bnd_label_new: label of the cur tp

% if cell_loc_map in detection region -> match
% else pick the nearest neighbor. if one of the cell < 0.8 -> division

movieInfo.tracks = {};
movieInfo.track_num = 0;
movieInfo.track_ids = zeros(8e6,1);
movieInfo.parents = cell(8e6,1);

im_sz = size(refine_res);
for tt = start_tp+1:end_tp
    tt
    cell_loc_map = mapCellLoc(cell_loc, tt-1, tt, driftInfo, drift);
    cell_size_old = cell_size;
    load(['/work/public/Embryo/Amat2014/Detection/Wei_refine_res/' sprintf('%05d', tt-1) '.mat']);
    [cell_loc, cell_size] = refineRes2cellLoc(refine_res);
    
    n_perframe(tt) = size(cell_loc,1);
    n_cumsum(tt+1) = n_cumsum(tt) + n_perframe(tt);
    frames(n_cumsum(tt)+1:n_cumsum(tt+1)) = tt;
    orgCoord(n_cumsum(tt)+1:n_cumsum(tt+1),:) = cell_loc .* [2 2 1];

    matching_j2i = nan(size(cell_loc,1),1);
    matching_i2j = nan(size(cell_loc_map,1),1);
    bnd_label_new = nan(size(cell_loc,1),1);
    idx = sub2ind_direct(im_sz, round(cell_loc_map(:,2)*2),...
                round(cell_loc_map(:,1)*2), round(cell_loc_map(:,3)));
    for ii = 1:length(idx)
        % try to match boundary cells as much as possible
        jj = refine_res(idx(ii));
        if jj > 0
            if isnan(matching_j2i(jj))
                matching_j2i(jj) = ii; 
                matching_i2j(ii) = jj;
                bnd_label_new(jj) = boundary_label(ii);

                parent_id = n_cumsum(tt-1)+ii;
                child_id = n_cumsum(tt)+jj;
                movieInfo = addChild(parent_id, child_id, movieInfo);

            elseif boundary_label(ii) == 1
                matching_j2i(jj) = ii;
                matching_i2j(ii) = jj;
                bnd_label_new(jj) = 1;

                parent_id = n_cumsum(tt-1)+ii;
                child_id = n_cumsum(tt)+jj;
                movieInfo = deleteChild(child_id, movieInfo);
                movieInfo = addChild(parent_id, child_id, movieInfo);
            end
        end
    end
    cell_dist = pdist2(cell_loc_map.*resolution, cell_loc.*resolution);
    for jj = 1:size(cell_loc,1)
        % further processing non matched cells
        if isnan(matching_j2i(jj))
            [~, cell_ord] = sort(cell_dist(:,jj));
            if isnan(matching_i2j(cell_ord(1))) && cell_dist(cell_ord(1), jj) < max_dist
                % case 1: nearest neighbor is not assigned
                matching_j2i(jj) = cell_ord(1);
                matching_i2j(cell_ord(1)) = jj;
                bnd_label_new(jj) = boundary_label(cell_ord(1));

                parent_id = n_cumsum(tt-1) + cell_ord(1);
                child_id = n_cumsum(tt)+jj;
                movieInfo = addChild(parent_id, child_id, movieInfo);
            elseif isnan(matching_i2j(cell_ord(2))) && cell_dist(cell_ord(2), jj) < max_dist
                % case 2: 2nd-nearest neighbor is not assigned
                matching_j2i(jj) = cell_ord(2);
                matching_i2j(cell_ord(2)) = jj;
                bnd_label_new(jj) = boundary_label(cell_ord(2));  

                parent_id = n_cumsum(tt-1) + cell_ord(2);
                child_id = n_cumsum(tt)+jj;
                movieInfo = addChild(parent_id, child_id, movieInfo);
            elseif cell_dist(cell_ord(1), jj) >= max_dist
                % not tracked
                bnd_label_new(jj) = 0;
            else
                % parent: cell_ord(1)
                % another_child = matching_i2j(cell_ord(1))
                size_ratio = min(cell_size(jj), cell_size(matching_i2j(cell_ord(1))))...
                            / cell_size_old(cell_ord(1));
                if size_ratio < size_ratio_thre
                    % case 3: matched neighbor or itself size ratio < 0.8
                    bnd_label_new(jj) = boundary_label(cell_ord(1));

                    parent_id = n_cumsum(tt-1) + cell_ord(1);
                    child_id = n_cumsum(tt)+jj;
                    movieInfo = addDivision(parent_id, child_id, movieInfo);
                elseif cell_dist(cell_ord(2), jj) >= max_dist
                    bnd_label_new(jj) = 0;
                else
                    % case 4: matched 2-nd neighbor or itself size ratio < 0.8
                    size_ratio = min(cell_size(jj), cell_size(matching_i2j(cell_ord(2))))...
                            / cell_size_old(cell_ord(2));
                    if size_ratio < size_ratio_thre
                        bnd_label_new(jj) = boundary_label(cell_ord(2));

                        parent_id = n_cumsum(tt-1) + cell_ord(2);
                        child_id = n_cumsum(tt)+jj;
                        movieInfo = addDivision(parent_id, child_id, movieInfo);
                    else
                        bnd_label_new(jj) = 0;
                    end
                end
            end
        end
    end
    boundary_label = bnd_label_new;
    cell_loc_all{tt} = [cell_loc boundary_label];
end
movieInfo.frames = movieInfo.frames(1:n_cumsum(end));
movieInfo.orgCoord = orgCoord(1:n_cumsum(end), :);
movieInfo.track_ids = frames(1:n_cumsum(end));
movieInfo.parents = movieInfo.parents(1:n_cumsum(end));
movieInfo.n_perframe = n_perframe;


function movieInfo = deleteChild(child_id, movieInfo)
    parent_id = movieInfo.parents{child_id};
    track_id = movieInfo.track_ids(child_id);
    if  movieInfo.track_ids(parent_id) ~= track_id
        error('Inconsist pair');
    end
    if movieInfo.tracks{track_id}(end) ~= child_id
        error('Last is not the child');
    end
    movieInfo.tracks{track_id}(end) = [];
    movieInfo.parents{child_id} = [];
    movieInfo.track_ids(child_id) = 0;
end

function movieInfo = addDivision(parent_id, child_id, movieInfo)
    % build a new track no matter exist or not
    movieInfo.track_num = movieInfo.track_num + 1;
    movieInfo.tracks{movieInfo.track_num} = [parent_id child_id];
    if movieInfo.track_ids(parent_id) == 0
        error('Division error');
    end
    movieInfo.track_ids(child_id) = movieInfo.track_num;
    movieInfo.parents{child_id} = parent_id;
end


function movieInfo = addChild(parent_id, child_id, movieInfo)
    if movieInfo.track_ids(parent_id) == 0
        movieInfo.track_num = movieInfo.track_num + 1;
        movieInfo.tracks{movieInfo.track_num} = [parent_id child_id];
        movieInfo.track_ids(parent_id) = movieInfo.track_num;
        movieInfo.track_ids(child_id) = movieInfo.track_num;
    else
        track_id = movieInfo.track_ids(parent_id);
        movieInfo.tracks{track_id}(end+1) = child_id;
        movieInfo.track_ids(child_id) = track_id;
    end
    movieInfo.parents{child_id} = parent_id;
end

function cell_loc_new = mapCellLoc(cell_loc, start_tp, end_tp, driftInfo, drift)
    cell_loc_new = cell_loc;
    for tt = start_tp:end_tp-1
        y_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.y_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        x_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.x_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        z_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.z_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        cell_loc_new = cell_loc_new + [y_bias x_bias z_bias];
    end
end

function [cell_loc, cell_size] = refineRes2cellLoc(refine_res)
    % find cell center from detection results
    prop = regionprops3(refine_res, 'VoxelList', 'VoxelIdxList');
    cell_size = cellfun(@length, prop.VoxelIdxList);
    cell_loc = cellfun(@mean, prop.VoxelList, 'UniformOutput', false);
    cell_loc = cell2mat(cell_loc);
    cell_loc(:,1) = cell_loc(:,1) / 2;
    cell_loc(:,2) = cell_loc(:,2) / 2;
end
