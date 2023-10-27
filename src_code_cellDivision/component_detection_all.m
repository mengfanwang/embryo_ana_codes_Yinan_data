clc;clear;clearvars -except movieInfo
% Another naive rule: check isolated regions in a detecion
%%
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');

%%
t = 192;
div_num = 0;
data_size = [960 960 181];
data_template = zeros(data_size);
cell_num = length(movieInfo.xCoord);
detect_splitPair = zeros(cell_num,2);
for cell_id = 1:cell_num
    if mod(cell_id,10000) == 0
        fprintf('%d / %d\n', cell_id, cell_num);     
    end
    [vBin, vIdx, ~, ~] = crop3D(data_template, movieInfo.voxIdx{cell_id}, [0 0 0]);
    vBin(ismember(vIdx, movieInfo.voxIdx{cell_id})) = 1;
    cc = bwconncomp(vBin, 6);
    if cc.NumObjects == 2 && ~isempty(movieInfo.parents{cell_id})
        size_ratio = cellfun(@length, cc.PixelIdxList);
        size_ratio = max(size_ratio) / min(size_ratio);
        if size_ratio < 3
            div_num = div_num+1;
            detect_splitPair(div_num,:) = [movieInfo.parents{cell_id} cell_id];
        end
    end
end
detect_splitPair = detect_splitPair(1:div_num,:);


