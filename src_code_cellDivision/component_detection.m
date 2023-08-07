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
        cc = bwconncomp(vBin, 6);
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

%%
function cell_match = findMatchCell(movieInfo, cell_idx, cell_time)
    cell_candi = find(movieInfo.frames == cell_time);
    cell_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(cell_candi)));
    if isempty(cell_match)
        cell_match = nan;
    else
        cell_match = cell_match + sum(movieInfo.n_perframe(1:cell_time-1));
    end
end