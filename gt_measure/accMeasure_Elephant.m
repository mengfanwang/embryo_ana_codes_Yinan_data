clc;clear;close all;

resolution = [1 1 1];
node_file = '/work/Nova/embryo_res_folder/mengfan_data_res/ground_truth/Elephant_detection_410_Spot.csv';
edge_file = '/work/Nova/embryo_res_folder/mengfan_data_res/ground_truth/MastodonTable-0_192-Link-new-Link.csv';
gt_file = '/ssd1/Mengfan/data/gt_train_all.csv';

% load ground truth
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
node = readmatrix(node_file);
node = node(3:end, [2 5 6 7 4]);
node = [node 30*ones(size(node,1),1)];
edge = readmatrix(edge_file);
edge = edge(3:end, [4 5]);
[tp, fp, fn, tp_list] = evaluation(gt_node, gt_edge, node, edge, resolution);
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));