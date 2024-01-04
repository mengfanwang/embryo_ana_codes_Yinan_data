clc;clear;close all;
load('/ssd1/Mengfan/result.mat');
node = [node 30*ones(size(node,1),1)]; % manual add radius
resolution = [1 1 5.86];

% get gt
% node: [id x y z t]
% edge: [id1 id2]
gt = readmatrix('/ssd1/Mengfan/data/gt_train.csv');
gt_node = gt(:,5:-1:1);
gt_node(:,2:4) = gt_node(:,2:4).*resolution;
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

%  function
addpath /home/mengfan/ForExecute/embryo_ana_codes_Yinan_data/CINDA/src_matlab

[tp, fp, fn] = evaluation(gt_node, gt_edge, node, edge, resolution);
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));