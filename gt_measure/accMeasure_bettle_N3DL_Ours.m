% clc;clear;close all;

% based on distance since the ground truth is sparse
% we use matlab coordiate, start from 1
t = 80;
gt_path = '/work/Mengfan/EmbryoData_other/Fluo-N3DL-TRIF_train/02_GT/TRA';

%% load ground truth
% get division
record = readmatrix(fullfile(gt_path, 'man_track.txt'));
record(:, 2:3) = record(:, 2:3) + 1;
% [id start_tp end_tp parent]
div_num = 0;
gt_division = zeros(size(record));
% [t parent child1 child2]
for ii = 1:size(record,1)
    if record(ii,4) > 0
        if ismember(record(ii,4), gt_division(:,2)) % already have one child
            parent_ind = find(record(ii,4) == gt_division(:,2));
            gt_division(parent_ind, 4) = record(ii,1);
        else % new division
            div_num = div_num + 1;
            parent_ind = record(ii,4);
            gt_division(div_num, 1) = record(parent_ind, 3);
            gt_division(div_num, 2) = parent_ind;
            gt_division(div_num, 3) = record(ii,1);
        end
    end
end
gt_division = gt_division(1:div_num, :);

% %%
% tif_files = dir(fullfile(gt_path, '*.tif'));
% tif_num = length(tif_files);
% gt_voxIdx = cell(tif_num, 1);
% max_track = 500;
% tic;
% for tt = 1:tif_num
%     tt
%     im = tifread(fullfile(gt_path, tif_files(tt).name));
%     im = imdilate(im, strel("sphere",2));
%     cell_ids = unique(im);
%     cell_ids = cell_ids(2:end);           
%     gt_voxIdx{tt} = cell(max_track,1);
%     for ii = 1:length(cell_ids)
% %         if broken_flag(cell_ids(ii)) == 0 
%         gt_voxIdx{tt}{cell_ids(ii)} = find(im == cell_ids(ii));
% %         end
%     end
% end
% toc

%%
load('bettle_N3DL_gtVox.mat'); max_track = 500; im_sz = [1820 1000 975];
n_perframe = zeros(80,1);
for tt = 1:t
    n_perframe(tt) = sum(~cellfun(@isempty, gt_voxIdx{tt}));
end
node_num = sum(n_perframe);
% cum_perframe = [0; cumsum(n_perframe)];
node_id = 0; edge_id = 0;
gt_node = zeros(node_num, 5); gt_edge = zeros(node_num*5, 2);
for tt = 1:t
%     node_id  = node_id + 1;
    for ii = 1:max_track
        if ~isempty(gt_voxIdx{tt}{ii})
            node_id  = node_id + 1;
            [xx, yy, zz] = ind2sub_direct(im_sz, round(mean(gt_voxIdx{tt}{ii})));
            gt_node(node_id, :) = [max_track*(tt-1)+ii yy xx zz tt];
            if tt > 1 && ~isempty(gt_voxIdx{tt-1}{ii})
                edge_id  = edge_id + 1;
                gt_edge(edge_id, :) = [max_track*(tt-2)+ii max_track*(tt-1)+ii];
            end
            div_id = find(ismember(record(:, [1 2]), [ii tt], 'rows'));
            if ~isempty(div_id) && (record(div_id, 4) > 0)
                edge_id = edge_id + 1;
                gt_edge(edge_id, :) = [max_track*(tt-2)+record(div_id, 4)  max_track*(tt-1)+ii];
            end
        end
    end
end
gt_edge = gt_edge(1:edge_id, :);

%% load detection result
% load('/work/Mengfan/EmbryoData_other/Fluo-N3DL-TRIF_train/02_Tracking/0208_000_079/movieInfo.mat');
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
resolution = [2 2 2];
[tp, fp, fn, tp_list] = evaluation(gt_node, gt_edge, node, edge, resolution);
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));
