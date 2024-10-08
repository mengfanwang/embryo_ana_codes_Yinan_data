clc;clear;close all;
% Another naive rule: check isolated regions in a detecion
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0504_000_191/movieInfo.mat');
%% get ground truth
addpath ../gt_measure/
path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(path, 1);
% link table: s t
% spot table: id t y x z
link_table = link_table(:,1:2);
spot_table = spot_table(:,1:5);

% convert to tracking result corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1/2 1/2 1/5.86] + 1;

%%
t = 192;
div_num = 0;
div_list = zeros(1000,1);
cc_list = cell(1000,1);
vid_list = cell(1000,2);
detect_splitPair = zeros(1000,3);
data_size = [960 960 181];
cell_num = size(spot_table,1);
for ii = 1:cell_num
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, cell_num);     
    end
    cell_loc = round(spot_table(ii, 3:5));
    cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
    cell_time = spot_table(ii,2);
    cell_match = findMatchCell(movieInfo, cell_idx, cell_time);
    if ~isnan(cell_match) && ~ismember(cell_match, detect_splitPair)
        [vidOut, vIdx, ~, ~] = crop3D(embryo_vid{cell_time}, movieInfo.voxIdx{cell_match}, [0 0 0]);
        vBin = zeros(size(vidOut));
        vBin(ismember(vIdx, movieInfo.voxIdx{cell_match})) = 1;
        cc = bwconncomp(vBin);
        if cc.NumObjects == 2 && ~isempty(movieInfo.parents{cell_match})
            size_ratio = cellfun(@length, cc.PixelIdxList);
            size_ratio = max(size_ratio) / min(size_ratio);
            if size_ratio < 3
                div_num = div_num+1;
                div_list(div_num) = ii;
                cc_list{div_num} = cc;
                vid_list{div_num,1} = vidOut;
                vid_list{div_num,2} = vBin;
                detect_splitPair(div_num,:) = [movieInfo.parents{cell_match} cell_match cell_match];
            end
        end
    end
end
div_list = div_list(1:div_num);
cc_list = cc_list(1:div_num);
detect_splitPair = detect_splitPair(1:div_num,:);
vid_list = vid_list(1:div_num,:);

%% get the location of cells 
a = spot_table(div_list,:);
a(:,2) = a(:,2) - 1;
a(:,3:5) = (a(:,3:5)-1) .* [2 2 1];
b = a;
b(:,5) = b(:,5)*5.86;

% 22-23 25-26 are FPs
c = cellfun(@(x) cellfun(@length, x.PixelIdxList), cc_list, 'UniformOutput', false)
c = cellfun(@(x) max(x(1),x(2))/min(x(1),x(2)), c)
%% check division detection accuracy
[TP, TN, FP, gt_divPair, detect_code_list] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) ;
fprintf('TP: %d  TN: %d  FP: %d\n', TP, TN, FP);
% 5: wrong association
% 11: FP detection
% 6 17 32 36 44: postpone 
% 26: sudden appear cell


%% convert FP to movieInfo format for visiualization
tp_flag = ismember(detect_divPair, gt_divPair, 'rows');
fp_pair = detect_divPair(~tp_flag,:);
fp_tracks = cell(size(fp_pair,1), 1);
movieInfo_temp = movieInfo;
movieInfo_temp.parents = cell(size(movieInfo.parents));
movieInfo_temp.kids = cell(size(movieInfo.kids));
for ii = 1:size(fp_pair,1)
    fp_tracks{ii} = fp_pair(ii,:);
    movieInfo_temp.parents{fp_pair(ii,2)} = fp_pair(ii,1);
    movieInfo_temp.kids{fp_pair(ii,1)} = fp_pair(ii,2);
end
movieInfo_temp.tracks = fp_tracks;
mat2tgmm(movieInfo_temp, '/work/Nova/embryo_res_folder/mengfan_data_res/temp');

%%

function [TP, TN, FP, gt_divPair, detect_code_list] = div_acc(detect_divPair, G, movieInfo, spot_table, data_size) 
    % get ground truth pairs
    d_out = outdegree(G);
    gt_list = G.Nodes.Name(find(d_out == 2));
    gt_divPair = nan(size(gt_list,1), 4);  %[parent_idx child_idx*2 parent_time]
    for ii = 1:length(gt_list)
        gt_pair = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_list{ii}),2);
        gt_pair = {gt_list{ii}; gt_pair{1}; gt_pair{2}};
        for jj = 1:3
            cell_ind = str2num(gt_pair{jj});
            cell_ind = find(spot_table(:,1) == cell_ind);
            cell_loc = round(spot_table(cell_ind, 3:5));
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            gt_divPair(ii,jj) = cell_idx;
            if jj == 1
                gt_divPair(ii, 4) = spot_table(cell_ind, 2);
            end
        end
    end
    % get postponed gt pairs
    gt_postPair = nan(size(gt_list,1), 4);
    for ii = 1:length(gt_list)
        gt_pair = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_list{ii}),2);
        for jj = 1:3
            if jj == 1
                cell_ind1 = str2num(gt_pair{1});
                cell_ind1 = find(spot_table(:,1) == cell_ind1);
                cell_loc1 = spot_table(cell_ind1, 3:5);
                cell_ind2 = str2num(gt_pair{2});
                cell_ind2 = find(spot_table(:,1) == cell_ind2);
                cell_loc2 = spot_table(cell_ind2, 3:5);
                cell_loc = round((cell_loc1 + cell_loc2)/2);
                gt_postPair(ii, 4) = gt_divPair(ii, 4)+1;
            else
                cell_ind = G.Edges.EndNodes(ismember(G.Edges.EndNodes(:,1), gt_pair{jj-1}),2);
                if length(cell_ind) == 1
                    cell_ind = find(spot_table(:,1) == str2num(cell_ind{1}));
                    cell_loc = round(spot_table(cell_ind, 3:5));
                elseif isempty(cell_ind)
                    cell_loc = nan(1,3);
                else
                    error('Consective division.');
                end
            end
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            gt_postPair(ii,jj) = cell_idx;
        end
    end
    % add missed annotation
    gt_addPair = [53044942 51210387 55802059 40;
                  116506827 113737225 117436109 60];
    gt_divPair = [gt_divPair; gt_addPair];
    gt_addPair(:,4) = gt_addPair(:,4) + 1;
    gt_postPair = [gt_postPair; gt_addPair];
    % remove one incorrect annotation
    error_loc = find(cellfun(@(x)strcmp(x, '34390'), gt_list));
    gt_list(error_loc) = [];
    gt_divPair(error_loc,:) = [];
    gt_postPair(error_loc,:) = [];
    gt_correctPair = [41064118 41989558 41059316 17;];
    gt_divPair = [gt_divPair; gt_correctPair];
    gt_correctPair(:,4) = gt_correctPair(:,4) + 1;
    gt_postPair = [gt_postPair; gt_correctPair];
  

    gt_code_list = zeros(size(gt_divPair,1), 1);
    detect_code_list = zeros(size(detect_divPair,1),1);
    % error code:
    % 0: missing parent
    % 1: correct
    % 2: child1 incorrect
    % 3: child2 incorrect
    % 4: both incorrect 
    for ii = 1:size(gt_divPair,1)
        parent_id = findMatchCell(movieInfo, gt_divPair(ii,1), gt_divPair(ii,4));
        detect_loc = find(parent_id == detect_divPair(:,1));
        if isnan(parent_id) || isempty(detect_loc)                          % case 0: missing parent
            error_code = 0;
        else
            child1_id = findMatchCell(movieInfo, gt_divPair(ii,2), gt_divPair(ii,4)+1);
            child2_id = findMatchCell(movieInfo, gt_divPair(ii,3), gt_divPair(ii,4)+1);
            if ismember(child1_id, detect_divPair(detect_loc,2:3)) && ...
               ismember(child2_id, detect_divPair(detect_loc,2:3))          % case 1: correct
                gt_code_list(ii) = 1;
                detect_code_list(detect_loc) = 1;
                continue;
            elseif ~ismember(child1_id, detect_divPair(detect_loc,2:3)) 
                error_code = 2;                                             % case 2: child1 incorrect
            elseif ~ismember(child2_id, detect_divPair(detect_loc,2:3)) 
                error_code = 3;                                             % case 3: child2 incorrect   
            else
                error_code = 4;                                             % case 4: both incorrect
            end
        end
        gt_code_list(ii) = error_code;
        if ~isempty(detect_loc)
            detect_code_list(detect_loc) = error_code;
        end
        % last try: postpone gt to 1 frame later
        parent_id = findMatchCell(movieInfo, gt_postPair(ii,1), gt_postPair(ii,4));
        detect_loc = find(parent_id == detect_divPair(:,1));
        child1_id = findMatchCell(movieInfo, gt_postPair(ii,2), gt_postPair(ii,4)+1);
        child2_id = findMatchCell(movieInfo, gt_postPair(ii,3), gt_postPair(ii,4)+1);
        if ismember(child1_id, detect_divPair(detect_loc,2:3)) && ...
           ismember(child2_id, detect_divPair(detect_loc,2:3))              % case 1: correct
            gt_code_list(ii) = 1;
            detect_code_list(detect_loc) = 1;
            continue;
        end
    end

    
    TP = sum(gt_code_list == 1);
    TN = length(gt_code_list) - TP;
    FP = size(detect_divPair,1) - TP;
end

function cell_match = findMatchCell(movieInfo, cell_idx, cell_time)
    cell_candi = find(movieInfo.frames == cell_time);
    cell_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(cell_candi)));
    if isempty(cell_match)
        cell_match = nan;
    else
        cell_match = cell_match + sum(movieInfo.n_perframe(1:cell_time-1));
    end
end