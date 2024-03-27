clc;clear;close all;

t = 80;
gt_path = '/work/Mengfan/EmbryoData_other/Fluo-N3DL-TRIF_train/02_GT';

%%
max_track = 500;
sc_f = [1 1 1];

node_num = max_track * t;
movieInfo.n_perframe = ones(t,1) * max_track;
movieInfo.xCoord = zeros(node_num,1);
movieInfo.yCoord = zeros(node_num,1);
movieInfo.zCoord = zeros(node_num,1);
movieInfo.frames = zeros(node_num,1);
movieInfo.tracks = cell(max_track,1);
movieInfo.parents = cell(node_num,1);

% get division
record = readmatrix(fullfile(gt_path, 'TRA', 'man_track.txt'));
record(:, 2:3) = record(:, 2:3) + 1;
% [id start_tp end_tp parent]
for ii = 1:size(record,1)
    if record(ii, 4) > 0
        tt = record(ii, 2);
        parent = record(ii,4) + (tt-2) * max_track;
        child = record(ii,1) + (tt-1) * max_track;
        movieInfo.tracks{record(ii,1)} = [parent child];
        movieInfo.parents{child} = parent;
    end
end

tif_files = dir(fullfile(gt_path, 'TRA', '*.tif'));
tic;
for tt = 1:t
    tt
    im = tifread(fullfile(gt_path, 'TRA', tif_files(tt).name));
    s = regionprops3(im, 'VoxelList');
    s = s.VoxelList;
    for ii = 1:length(s)
        if ~isempty(s{ii})
            id = ii + (tt-1) * max_track;
            movieInfo.frames(id) = tt;
            movieInfo.xCoord(id) = mean(s{ii}(:,1));
            movieInfo.yCoord(id) = mean(s{ii}(:,2));
            movieInfo.zCoord(id) = mean(s{ii}(:,3));
            if isempty(movieInfo.parents{id})
                if tt > 1
                    parent_id = ii + (tt-2) * max_track;
                    movieInfo.parents{id} = parent_id;
                end
                movieInfo.tracks{ii} = [movieInfo.tracks{ii} id];
            end
        end
    end
end
toc

mat2tgmm(movieInfo, fullfile(gt_path, 'tgmm_format'), sc_f);