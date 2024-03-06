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
table(:, 4) = table(:,4) .* 5;

t = end_tp - start_tp + 1;
node_num = size(table, 1);
edge_num = 0;
gt_edge = zeros(node_num,2);
for ii = 1:node_num
    parent_id = find(table(ii,5) == table(:,1));
    if ~isempty(parent_id)
        edge_num = edge_num + 1;
        gt_edge(edge_num, :) = [parent_id ii];
    end
end
gt_edge = gt_edge(1:edge_num, :);
gt_node = [(1:node_num)' table(:, [3 2 4 6])];

%%
% load detection data
% node:[id x y z t radius]
num_node = length(movieInfo.xCoord);
node = zeros(num_node, 6);
node(:,1) = (1:num_node)';
node(:,2) = movieInfo.xCoord;
node(:,3) = movieInfo.yCoord;
node(:,4) = movieInfo.zCoord;
node(:,5) = movieInfo.frames;
node(:,6) = 30; 
edge = zeros(num_node, 2);
edge_num = 0;
for ii = 1:length(movieInfo.tracks)
    for jj = 1:length(movieInfo.tracks{ii})-1
        edge_num = edge_num + 1;
        edge(edge_num,:) = [movieInfo.tracks{ii}(jj) movieInfo.tracks{ii}(jj+1)]; 
    end
end
edge = edge(1:edge_num, :);

%%
resolution = [2 2 5];
[tp, fp, fn, tp_list] = evaluation(gt_node, gt_edge, node, edge, resolution);
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));
