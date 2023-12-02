clc;clear;close all;
% merge broken tracks as much as possible for lineage reconstruction
% load('/work/Nova/embryo_res_folder/mengfan_data_res/merge_division/movieInfo.mat');
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');
movieInfo_temp = movieInfo;

%% initial g
g = graphPara_cell(length(movieInfo.xCoord));
g.loopCnt = 100;

g.translation_path = '/work/Mengfan/Embryo/Registration/nonrigid_result_view10_l5_s1_linear';
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
g.driftInfo.z_batch = z_batch; 

%% find neighbor
dbstop if error
% parameter setting
paras.im_resolution = [1 1 5.86];
paras.max_dist = 50;
paras.max_nei = 5;

% find track head and tails
% track_heads = find(cellfun(@length, movieInfo.parents)==0);
% track_heads(track_heads>length(movieInfo.xCoord)) = [];
% track_tails = find(cellfun(@length, movieInfo.kids)==0);
% track_tails(track_tails>length(movieInfo.xCoord)) = [];
% valid_ids = union(track_heads, track_tails);
curTracks = movieInfo.tracks;
curTracks(cellfun(@length, curTracks)==0) = [];
track_heads = []; track_tails = [];
for ii = 1:length(curTracks)
    track_heads(end+1) = curTracks{ii}(1);
    track_tails(end+1) = curTracks{ii}(end);
end

% for ii = 1:length(valid_ids)
%     if mod(ii, 1000) == 0
%         fprintf('%d / %d\n', ii, length(valid_ids));
%     end
%     if movieInfo.frames(valid_ids(ii)) < length(movieInfo.n_perframe)
%         movieInfo.nei{valid_ids(ii)} = findNeighbor(movieInfo, valid_ids(ii), paras)';
%     end
%     if movieInfo.frames(valid_ids(ii)) > 1
%         movieInfo.preNei{valid_ids(ii)} = findNeighborPrev(movieInfo, valid_ids(ii), paras)';
%     end
% end
 
mergeTracks = curTracks;
for ii = 1:length(track_tails)
    if mod(ii, 1000) == 0
        fprintf('%d / %d\n', ii, length(track_tails));
    end
    if movieInfo.frames(track_tails(ii)) < length(movieInfo.n_perframe)
        candidate = findNeighbor(movieInfo, track_tails(ii), paras);
        candidate = candidate(ismember(candidate, track_heads));
        if ~isempty(candidate)
            tail = find(track_heads == candidate(1));
            mergeTracks{ii} = [mergeTracks{ii}; mergeTracks{tail(1)}];
            mergeTracks{tail(1)} = [];
        end
    end
end
mergeTracks(cellfun(@length, mergeTracks)==0) = [];

%% merge broken tracks
dbstop if error
movieInfo = mergeBrokenTracksV3(movieInfo, g);

function [neighbors, dist2nei] = findNeighbor(movieInfo, cur_i, paras)
    % find the nearest neighbors in the next frame with given conditions  
    im_resolution = paras.im_resolution;
    max_dist = paras.max_dist;
    max_nei = paras.max_nei;
    
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
    dist2nei = neighbor_dist(neighbor_dist < max_dist);
end

function [neighbors, dist2nei] = findNeighborPrev(movieInfo, cur_i, paras)
    % find the nearest neighbors in the next frame with given conditions  
    im_resolution = paras.im_resolution;
    max_dist = paras.max_dist;
    max_nei = paras.max_nei;
    
    tt = movieInfo.frames(cur_i);
    curCentroid = movieInfo.orgCoord(cur_i,:);
    prevCentroid = movieInfo.orgCoord(movieInfo.frames==(tt-1),:);
    
    % candidate neighbors roughly selection
    drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt-1, tt, movieInfo.drift, movieInfo.driftInfo);
    dist_pair = pdist2(curCentroid.*im_resolution, (prevCentroid + drift).*im_resolution); % not accurate but enough
    [neighbor_dist, neighbor_candidate] = sort(dist_pair);
    neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
    neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
    neighbors = neighbor_candidate(neighbor_dist < max_dist);
    neighbors = neighbors + sum(movieInfo.n_perframe(1:tt-2));
    dist2nei = neighbor_dist(neighbor_dist < max_dist);
end