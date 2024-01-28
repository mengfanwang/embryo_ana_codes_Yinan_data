clc;clear;close all;

%% load ground truth file
gt_folder = '/work/Mengfan/Embryo/Amat2014/ground_truth';
start_tp = 300; end_tp = 399;

gt_files = dir(fullfile(gt_folder, '*.mat'));
table = [];
for ii = 1:length(gt_files)
    % 'node id', 'type', 'x', 'y', 'z', 'radius', 'parent node id', 'timepoint', 'tag', 'lineage id'
    load(fullfile(gt_folder, gt_files(ii).name));
    tracks = tracks(tracks(:,8) <= end_tp & tracks(:,8) >= start_tp, :);
    tracks(:, 3:4) = tracks(:, 3:4)/333;
    tracks(:, 5) = tracks(:, 5)/2985;
    table = [table; tracks];
end
% 'node id', 'x', 'y', 'z', 'parent node id', 'timepoint'
table = table(:, [1 3 4 5 7 8]);
table = sortrows(table, 6);
table(:, 6) = table(:,6) - start_tp + 1;
table(:, 2:4) = table(:,2:4) + 1;   % for matlab

t = end_tp - start_tp + 1;
node_num = size(table, 1);
n_perframe = zeros(t,1);
for tt = 1:t
    n_perframe(tt) = sum(table(:,6) == tt);
end
movieInfo.xCoord = table(:,3);
movieInfo.yCoord = table(:,2);
movieInfo.zCoord = table(:,4);
tracks = {};
parents = cell(node_num,1);
kids = cell(node_num,1);
track_ids = zeros(node_num, 1);
track_num = 0;
for ii = 1:node_num
    parent_id = find(table(ii,5) == table(:,1));
    if isempty(parent_id)
        track_num = track_num + 1;
        tracks{track_num} = ii;
        track_ids(ii) = track_num;
    else
        if isempty(kids{parent_id})
            track_id = track_ids(parent_id);
            tracks{track_id} = [tracks{track_id}; ii];
            track_ids(ii) = track_id;
            parents{ii} = parent_id;
            kids{parent_id} = ii;
        else
            % division
            track_num = track_num + 1;
            tracks{track_num} = [parent_id; ii];
            track_ids(ii) = track_num;
            parents{ii} = parent_id;
            kids{parent_id} = [kids{parent_id}; ii];
        end
    end
end
movieInfo.tracks = tracks;
movieInfo.parents = parents;
movieInfo.n_perframe = n_perframe;