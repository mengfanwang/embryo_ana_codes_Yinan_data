% given two batch of cell segmentation results, combine the mat files into
% unifiled mat files.

batch_path = cell(2, 1);
batch_path{1} = '/work/public/sameViewFusion/sameViewDetection_230-239_08';
batch_path{2} = '/work/public/sameViewFusion/sameViewDetection_240-249_08';

write_folder = '/work/public/sameViewFusion/TwentyFrames/sameViewDetection_230-249_08';

for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d');
    else
        cur = load(fullfile(batch_path{i}, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d');
        eig_res_2d = cat(1, eig_res_2d, cur.eig_res_2d);
        eig_res_3d = cat(1, eig_res_3d, cur.eig_res_3d);
    end
end
eig_overlay = cell(numel(eig_res_3d), 1);
save(fullfile(write_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d', 'eig_overlay', '-v7.3');
%clear eig_res_2d eig_res_3d eig_overlay


for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'synQuant_refine_res.mat'));
    else
        cur = load(fullfile(batch_path{i}, 'synQuant_refine_res_4d_v9.mat'));
        refine_res = cat(1, refine_res, cur.refine_res);
        threshold_res = cat(1, threshold_res, cur.threshold_res);
    end
end
save(fullfile(write_folder, 'synQuant_refine_res_4d_v9.mat'), 'refine_res', 'threshold_res','-v7.3');
%clear refine_res threshold_res

for i = 1:numel(batch_path)
    if i == 1
        load(fullfile(batch_path{i}, 'varianceMap.mat'));
    else
        cur = load(fullfile(batch_path{i}, 'varianceMap.mat'));
        varMap = cat(1, varMap, cur.varMap);
    end
end
save(fullfile(write_folder, 'varianceMap.mat'), 'varMap','-v7.3');
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

save(fullfile(write_folder, 'synQuant_res.mat'), 'fMaps', 'id_mat', 'z_mat', '-v7.3');
%clear fMaps id_mat z_mat cur


