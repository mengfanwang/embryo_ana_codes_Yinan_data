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


gt_detection_folder = fullfile(gt_path, 'track_identities/TIF');
tif_files = dir(fullfile(gt_detection_folder, '*.tiff'));
tif_num = length(tif_files);
im = cell(tif_num,1);
cell_ids = cell(tif_num,1);
n_perframe = zeros(tif_num,1);
n_cumsum = zeros(tif_num+1,1);
orgCoord = zeros(1000000,3);
frames = zeros(1000000,1);
max_track = 0;
for tt = 1:tif_num
    tt
    im{tt} = tifread(fullfile(gt_detection_folder, tif_files(tt).name));
    im{tt}(im{tt} == 65535) = 0;
    im_sz = size(im{tt});
    cell_ids{tt} = unique(im{tt});
    cell_ids{tt} = cell_ids{tt}(2:end);
    n_perframe(tt) = length(cell_ids{tt});
    max_track = max(cell_ids{tt});
    n_cumsum(tt+1) = sum(n_perframe(1:tt));
    for ii = 1:n_perframe(tt)
        [yy, xx, zz] = ind2sub_direct(im_sz, find(im{tt} == cell_ids{tt}(ii)));
        orgCoord(n_cumsum(tt) + ii, :) = mean([xx yy zz]);
    end
    frames(n_cumsum(tt)+1:n_cumsum(tt+1)) = tt;
end
cell_num = n_cumsum(end);
frames = frames(1:cell_num);
orgCoord = orgCoord(1:cell_num,:);

%% build movieInfo file
movieInfo.xCoord = orgCoord(:,1);
movieInfo.yCoord = orgCoord(:,2);
movieInfo.zCoord = orgCoord(:,3);
movieInfo.n_perframe = n_perframe;
movieInfo.parents = cell(cell_num, 1);
movieInfo.tracks = cell(max_track, 1);
broken_flag = zeros(max_track, 1);
for tt = 1:tif_num - 1
    for ii = 1:n_perframe(tt)
        track_id = cell_ids{tt}(ii);
        % check division
        if ismember(track_id, gt_division(:,5)) || ismember(track_id, gt_division(:,7)) 
            if ismember(track_id, gt_division(:,5))
                parent_id = gt_division(track_id == gt_division(:,5),3);
            elseif ismember(track_id, gt_division(:,7))
                parent_id = gt_division(track_id == gt_division(:,7),3);
            end
            if ismember(parent_id, cell_ids{tt-1})
                movieInfo.tracks{track_id} = [n_cumsum(tt-1) + find(parent_id == cell_ids{tt-1}) ...
                                            n_cumsum(tt) + ii];
                movieInfo.parents{n_cumsum(tt) + ii} = n_cumsum(tt-1) + n_cumsum(tt-1) + find(parent_id == cell_ids{tt-1}); 
            end
       end
        if ismember(track_id, cell_ids{tt+1})
            if ~isempty(movieInfo.tracks{track_id})
                if frames(movieInfo.tracks{track_id}(end)) ~= tt
                    broken_flag(track_id) = 1;
                    fprintf('Error! Disconnect tracks. Count = %d\n', sum(broken_flag));               
                end
            else
                movieInfo.tracks{track_id} = n_cumsum(tt) + ii;
            end
            movieInfo.tracks{track_id}(end+1) = n_cumsum(tt+1) + find(track_id == cell_ids{tt+1});
            movieInfo.parents{n_cumsum(tt+1) + find(track_id == cell_ids{tt+1})} = n_cumsum(tt) + ii;
        end
    end
end
