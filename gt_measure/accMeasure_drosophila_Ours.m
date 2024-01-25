clc;clear;close all;

% please note that one cell may have multiple ids in the drosophila data
% we use matlab coordinate, start from 1
resolution = [1 1 1];
gt_path = '/work/Mengfan/EmbryoData_other/drosophila-cell-tracking/manual_annotations/convenience_files';

%% load ground truth
gt_division_file = fullfile(gt_path, 'CSV/manual_annotations_divisions.csv');
% [t parent_id/track child1/2 id_track]
gt_division = readmatrix(gt_division_file);
gt_division(:,1) = gt_division(:,1) + 1;

load('drosophila_broken_tracks.mat');
gt_detection_folder = fullfile(gt_path, 'track_identities/TIF');
tif_files = dir(fullfile(gt_detection_folder, '*.tiff'));
tif_num = length(tif_files);
gt_voxIdx = cell(tif_num, 1);
max_track = 3000;
tic;
for tt = 1:tif_num
    tt
    im = tifread(fullfile(gt_detection_folder, tif_files(tt).name));
    im(im == 65535) = 0;
    im = imdilate(im, strel("sphere",2));
    cell_ids = unique(im);
    cell_ids = cell_ids(2:end); 
    gt_voxIdx{tt} = cell(max_track,1);
    for ii = 1:length(cell_ids)
        if broken_flag(cell_ids(ii)) == 0 
            gt_voxIdx{tt}{cell_ids(ii)} = find(im == cell_ids(ii));
        end
    end
end
toc

%% check link num
parent_ids = cell(tif_num-1, 1);
for tt = 1:tif_num-1
    parent_ids{tt} = zeros(length(gt_voxIdx{tt}), 1);
    parent_cnt = 0;
    for ii = 1:length(gt_voxIdx{tt})
        if ~isempty(gt_voxIdx{tt}{ii}) && ~isempty(gt_voxIdx{tt+1}{ii})
            parent_cnt = parent_cnt + 1;
            parent_ids{tt}(parent_cnt) = ii;
        end
    end
    parent_ids{tt} = parent_ids{tt}(1:parent_cnt);
end
gt_edge_num = sum(cellfun(@length, parent_ids)) + size(gt_division,1);

%% load detection result
load('/work/Mengfan/EmbryoData_other/drosophila-cell-tracking/Our/Tracking/0103_0_99/movieInfo.mat');
% node:[id centerIdx t]
node_num = length(movieInfo.xCoord);
node = zeros(node_num, 3);
node(:,1) = (1:node_num)';
node(:,2) = sub2ind_direct(size(im), round(movieInfo.yCoord), ...
                    round(movieInfo.xCoord), round(movieInfo.zCoord));
node(:,3) = movieInfo.frames;
edge = zeros(node_num, 2);
edge_num = 0;
for ii = 1:length(movieInfo.tracks)
    for jj = 1:length(movieInfo.tracks{ii})-1
        edge_num = edge_num + 1;
        edge(edge_num,:) = [movieInfo.tracks{ii}(jj) movieInfo.tracks{ii}(jj+1)]; 
    end
end
edge = edge(1:edge_num, :);

%% evaluation
tp = 0; fp = 0; 
tp_division = 0;
tp_list = zeros(edge_num,1); fp_list = zeros(edge_num, 1);
detect2gt = zeros(node_num,1);
for ii = 1:edge_num
    if mod(ii, 1000) == 0
        fprintf('%d/%d\n', ii, edge_num);
    end
    tt = node(edge(ii,1),3);
    parent_id = 0;
    child_id = 0;
    for jj = 1:length(gt_voxIdx{tt})
        if ismember(node(edge(ii,1),2), gt_voxIdx{tt}{jj})
            parent_id = jj;
        end
    end
    for jj = 1:length(gt_voxIdx{tt+1})
        if ismember(node(edge(ii,2),2), gt_voxIdx{tt+1}{jj})
            child_id = jj;
        end
    end

    detect2gt(edge(ii,1)) = parent_id;
    detect2gt(edge(ii,2)) = child_id;

    if parent_id == 0
        if child_id > 0 && ~isempty(gt_voxIdx{tt}(child_id))
            fp = fp + 1;
            fp_list(ii) = 1;
        end
    else
        if child_id == parent_id
            tp = tp + 1;
            tp_list(ii) = 1;
        else
            % check division
            div_flag = 0;
            for jj = 1:size(gt_division,1)
                if gt_division(jj,1) == tt && gt_division(jj,3) == parent_id
                    if gt_division(jj,5) == child_id || gt_division(jj,7) == child_id
                        div_flag = 1;
                        break;
                    end
                end
            end
            if div_flag
                tp = tp + 1;
                tp_list(ii) = 1;
                tp_division = tp_division + 1;
            elseif ~isempty(gt_voxIdx{tt+1}(parent_id))
                fp = fp + 1;
                fp_list(ii) = 1;
            end
        end
    end
end
a = [detect2gt(edge) movieInfo.frames(edge(:,1)) movieInfo.frames(edge(:,2))];
[b, ib] = unique(a, 'rows');
tp = sum(tp_list(ib));
fp = sum(fp_list(ib));
fn = gt_edge_num - tp;
fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));