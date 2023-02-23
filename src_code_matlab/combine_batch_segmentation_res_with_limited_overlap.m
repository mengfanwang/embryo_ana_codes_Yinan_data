% given two batch of cell segmentation results, combine the mat files into
% two new batches. These two batches has |overlap| frames overlap.

batch_path = cell(2, 1);
% batch_path{1} = '/work/public/sameViewFusion/sameViewDetection_230-239_08';
% batch_path{2} = '/work/public/sameViewFusion/sameViewDetection_240-249_08';
batch_path{1} = '/work/public/sameViewFusion/sameViewDetection_050-069_11';
batch_path{2} = '/work/public/sameViewFusion/sameViewDetection_070-088_11';

overlap = 5; % e.g. combine timepoints 50-74 and 70-89 with 5 frames overlapped
write_folder{1} = '/work/public/sameViewFusion/merged_batch/sameViewDetection_050-074_11';
if ~exist(write_folder{1}) mkdir(write_folder{1}); end
% !!! -> NOTE: this is wrong! This script indeed generate 65-88 totally 24 time
% points
write_folder{2} = '/work/public/sameViewFusion/merged_batch/sameViewDetection_066-088_11';
if ~exist(write_folder{2}) mkdir(write_folder{2}); end

%% start combination
batch1_end = 25;
batch2_start = 21;
for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d');
        batch1_end = numel(eig_res_2d) + overlap; % 25
        batch2_start = numel(eig_res_2d) - overlap + 1; % 16 ==> time point from 65->88
    else
        cur = load(fullfile(batch_path{i}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d');
        eig_res_2d = cat(1, eig_res_2d, cur.eig_res_2d);
        eig_res_3d = cat(1, eig_res_3d, cur.eig_res_3d);
    end
end
eig_res_2d_all = eig_res_2d;
eig_res_3d_all = eig_res_3d;

eig_res_2d = eig_res_2d_all(1:batch1_end);
eig_res_3d = eig_res_3d_all(1:batch1_end);
eig_overlay = cell(numel(eig_res_3d), 1);
save(fullfile(write_folder{1}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d', 'eig_overlay', '-v7.3');
eig_res_2d = eig_res_2d_all(batch2_start:end);
eig_res_3d = eig_res_3d_all(batch2_start:end);
eig_overlay = cell(numel(eig_res_3d), 1);
save(fullfile(write_folder{2}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d', 'eig_overlay', '-v7.3');
%clear eig_res_2d eig_res_3d eig_overlay

%% less time consuming
for i = 1:numel(batch_path)
    if i == 1
        if exist(fullfile(batch_path{i}, 'synQuant_refine_res_4d_v9.mat'))
            load(fullfile(batch_path{i}, 'synQuant_refine_res_4d_v9.mat'));
        else
            load(fullfile(batch_path{i}, 'synQuant_refine_res.mat'));
        end
    else
        if exist(fullfile(batch_path{i}, 'synQuant_refine_res_4d_v9.mat'))
            cur = load(fullfile(batch_path{i}, 'synQuant_refine_res_4d_v9.mat'));
        else
            cur = load(fullfile(batch_path{i}, 'synQuant_refine_res.mat'));
        end
        refine_res = cat(1, refine_res, cur.refine_res);
        threshold_res = cat(1, threshold_res, cur.threshold_res);
    end
end
refine_res_all = refine_res;
threshold_res_all = threshold_res;

refine_res = refine_res_all(1:batch1_end);
threshold_res = threshold_res_all(1:batch1_end);
save(fullfile(write_folder{1}, 'synQuant_refine_res.mat'), 'refine_res', 'threshold_res','-v7.3');
refine_res = refine_res_all(batch2_start:end);
threshold_res = threshold_res_all(batch2_start:end);
save(fullfile(write_folder{2}, 'synQuant_refine_res.mat'), 'refine_res', 'threshold_res','-v7.3');
%clear refine_res threshold_res

for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'varianceMap.mat'));
    else
        cur = load(fullfile(batch_path{i}, 'varianceMap.mat'));
        varMap = cat(1, varMap, cur.varMap);
    end
end
varMap_all = varMap;
varMap = varMap_all(1:batch1_end);
save(fullfile(write_folder{1}, 'varianceMap.mat'), 'varMap','-v7.3');
varMap = varMap_all(batch2_start:end);
save(fullfile(write_folder{2}, 'varianceMap.mat'), 'varMap','-v7.3');
%clear varMap


for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'synQuant_res.mat'));
    else
        cur = load(fullfile(batch_path{i}, 'synQuant_res.mat'));
        fMaps = cat(1, fMaps, cur.fMaps);
        id_mat = cat(1, id_mat, cur.id_mat);
        z_mat = cat(1, z_mat, cur.z_mat);
    end
end
fMaps_all = fMaps;
id_mat_all = id_mat;
z_mat_all = z_mat;

fMaps = fMaps_all(1:batch1_end);
id_mat = id_mat_all(1:batch1_end);
z_mat = z_mat_all(1:batch1_end);
save(fullfile(write_folder{1}, 'synQuant_res.mat'), 'fMaps', 'id_mat', 'z_mat', '-v7.3');
fMaps = fMaps_all(batch2_start:end);
id_mat = id_mat_all(batch2_start:end);
z_mat = z_mat_all(batch2_start:end);
save(fullfile(write_folder{2}, 'synQuant_res.mat'), 'fMaps', 'id_mat', 'z_mat', '-v7.3');
%clear fMaps id_mat z_mat cur


