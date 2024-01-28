clc;clear;close all;

resolution = [2 2 5.86];
tic;
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');
toc
gt_file = '/ssd1/Mengfan/data/gt_train_all.csv';

%% load ground truth
gt = readmatrix(gt_file);
gt_node = gt(:,5:-1:1);
gt_edge = zeros(size(gt,1), 2);
edge_num = 0;
for ii = 1:size(gt,1)
    cell_id = gt(ii,5);
    parent_id = gt(ii,6);
    if ismember(parent_id, gt(:,5))
        edge_num = edge_num + 1;
        gt_edge(edge_num, :) = [parent_id cell_id];
    end
end
gt_edge = gt_edge(1:edge_num, :);

% load detection data
% node:[id x y z t radius]
num_node = length(movieInfo.xCoord);
node = zeros(num_node, 6);
node(:,1) = (1:num_node)';
node(:,2) = movieInfo.xCoord-1;
node(:,3) = movieInfo.yCoord-1;
node(:,4) = movieInfo.zCoord+1;
node(:,5) = movieInfo.frames-1;
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

% load division edges
load('/work/Nova/embryo_res_folder/mengfan_data_res/merge_division/detect_divPair.mat');
division_edge = [];
for ii = 1:size(detect_divPair,1)
    if ~ismember(detect_divPair(ii,[1 2]), edge, 'rows')
        division_edge = [division_edge; detect_divPair(ii,[1 2])];
    end
    if ~ismember(detect_divPair(ii,[1 3]), edge, 'rows')
        division_edge = [division_edge; detect_divPair(ii,[1 3])];
    end
end
edge = [edge; division_edge];

[tp, fp, fn, tp_list] = evaluation(gt_node, gt_edge, node, edge, resolution);
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));
