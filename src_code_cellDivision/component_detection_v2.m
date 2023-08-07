clc;clear;close all;
% Another naive rule: check isolated regions in a detecion
% v2: Not only consider isolated regions but also conncect cells
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

%% result
% 1-21: 0 eigenvalue
% 22: over-merge, not division
% 23: over-split, the cell has a strong gap
% 24: division
% 25: over-merge

% consider intensity total 100 completely seprated cells
% 1-31 are seperated   7 new found
% 3 5 7 14 31: division
% 29: over-merge, not division
% 32: division
% 33 34: very like to be division. not annotated, not sure
% 35: division
% 36: over-merge
% 37: division
% 38: over-merge. One part is the child of another division
% 39 40: division
% 41: false positive. Divide a cell into two parts. (Reserve high inten?)
% 42 43 44: division
% 45: over-merge. The child is wrongly seged with a part of another cell
% 46: very like to be division. not annotated, not sure
% 47: over-merge. One part is the child of another division
% 48: flase postivie. Very noise region.
% 49: false positive. Divide a cell into two parts. (Reserve high inten?)
% 50: over-merge


%%
addpath /home/mengfan/ForExecute/embryo_ana_codes_Yinan_data/src_code_cellSegment
load("cell_match_list.mat");
t = 192;
cell_num = size(spot_table,1);
vid_list = cell(cell_num,2);
detect_splitPair_all = nan(cell_num,3);
data_size = [960 960 181];
eigval_list = nan(cell_num,1);

tic;
for ii = 1:cell_num
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, cell_num);     
    end

    %
    cell_loc = round(spot_table(ii, 3:5));
    cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
    cell_time = spot_table(ii,2);
%     cell_match = findMatchCell(movieInfo, cell_idx, cell_time);
    cell_match = cell_match_list(ii);
    if ~isnan(cell_match) && ~ismember(cell_match, detect_splitPair_all)
        [vidOut, vIdx, ~, ~] = crop3D(embryo_vid{cell_time}, movieInfo.voxIdx{cell_match}, [0 0 0]);
        vBin = zeros(size(vidOut));
        vBin(ismember(vIdx, movieInfo.voxIdx{cell_match})) = 1;

        % using fiedler value and vector
        % build graph
        v_ind = find(vBin);
        node_reverse = zeros(max(v_ind),1);
        node_num = length(v_ind);
        edges = [];
        NodeNames = cell(node_num,1);
        XData  = zeros(node_num,1);
        YData  = zeros(node_num,1);
        ZData  = zeros(node_num,1);
        for jj = 1:node_num
            [~, ~, nei_ind] = neighbours(v_ind(jj), vBin, 6);
            nei_ind = nei_ind(ismember(nei_ind, v_ind));
            edges = [edges; [repmat(v_ind(jj), length(nei_ind), 1) nei_ind]];
            NodeNames{jj} = num2str(v_ind(jj));
            [XData(jj), YData(jj), ZData(jj)] = ind2sub_direct(size(vBin), v_ind(jj));
            node_reverse(v_ind(jj)) = jj;
        end
        edges = sort(edges,2);
        edges = unique(edges, 'rows');
        edge_num = size(edges,1);
        EndNodes = cell(edge_num,2);
        Weight = zeros(edge_num,1);
        for jj = 1:edge_num
            EndNodes{jj,1} = num2str(edges(jj,1));
            EndNodes{jj,2} = num2str(edges(jj,2));
            Weight(jj) = (vidOut(edges(jj,1)) + vidOut(edges(jj,2)))/2 + 1;
        end
        EdgeTable = table(EndNodes, Weight);
        NodeTable = table(NodeNames,'VariableNames',{'Name'});
        G = graph(EdgeTable, NodeTable);

%         L = laplacian(G);
        W = zeros(node_num, node_num);
        for jj = 1:edge_num
            W(node_reverse(edges(jj,1)),node_reverse(edges(jj,2))) = Weight(jj);
            W(node_reverse(edges(jj,2)),node_reverse(edges(jj,1))) = Weight(jj);
        end
        D = zeros(node_num, node_num);
        D_sqrtinv = zeros(node_num, node_num);
        for jj = 1:node_num
            D(jj,jj) = sum(W(jj,:));
            D_sqrtinv(jj,jj) = 1/sqrt(D(jj,jj));
        end
        L = D_sqrtinv*(D-W)*D_sqrtinv;
        [L_eigvec,L_eigval] = eigs(L,2,'smallestreal');
        % normalization to balance the size difference
        L_eigvec = D_sqrtinv*L_eigvec(:,2);
        L_eigval = L_eigval(2,2) / norm(L_eigvec) * sqrt(node_num);
        
        vec_sign = logical(kmeans(L_eigvec, 2) - 1);
        
%         % plot graph
%         node_color = vec_sign*[0 0.4470 0.7410] + (1-vec_sign)*[0.8500 0.3250 0.0980];
%         plot(G,'XData', XData, 'YData', YData, 'ZData', ZData,...
%             'NodeColor',node_color,'MarkerSize',5);

        vBin_pos = zeros(size(vBin));
        vBin_pos(v_ind(vec_sign)) = 1;
        cc_pos = bwconncomp(vBin_pos);
        vBin_neg = zeros(size(vBin));
        vBin_neg(v_ind(~vec_sign)) = 1;
        cc_neg = bwconncomp(vBin_neg);       

        if cc_pos.NumObjects == 1 && cc_neg.NumObjects == 1 ...
                && ~isempty(movieInfo.parents{cell_match})
            size_ratio = [length(cc_pos.PixelIdxList{1}) length(cc_neg.PixelIdxList{1})];
            size_ratio = max(size_ratio) / min(size_ratio);
            if size_ratio < 3
                vid_list{ii,1} = vidOut;
                vid_list{ii,2} = vBin;
                detect_splitPair_all(ii,:) = [movieInfo.parents{cell_match} cell_match cell_match];

                eigval_list(ii) = L_eigval;
            end
        end
    end
end 
toc

%% 3-sigma method, keep small eigen val cases
sep_flag = eigval_list < 1e-8;
small_eig_flag = ~sep_flag & eigval_list <...
    mean(eigval_list, 'omitnan') - 3*std(eigval_list, 'omitnan');
detect_splitPair{1} = detect_splitPair_all(sep_flag,:);
detect_splitPair{2} = detect_splitPair_all(small_eig_flag,:);


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