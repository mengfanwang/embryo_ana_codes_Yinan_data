% crop the original image data and also the cell segmentation results
% if intermediate results like principal curvature, variance exist, also 
% crop them.
home_folder = fullfile('/home', getenv('USER'));
addpath(fullfile(home_folder, 'Dropbox/cc_ImHandle/'));
save_folder = fullfile(home_folder, 'Desktop/embryo_res_folder');

image_folder_path = '/work/public/sameViewFusion/samViewFusion_240-249_08';
seg_res_folder_path = '/work/public/sameViewFusion/sameViewDetection_240-249_08';

% write_folder = '/home/ccwang/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/Mengfan_data_10012022/cropped';
% scale = [477, 1643; 247, 1313; 131, 170];

write_folder = '/home/ccwang/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/Mengfan_data_10012022/cropped_further';
scale = [560, 1059; 345, 844; 131, 170];
%% crop images
tif_files = dir(fullfile(image_folder_path, '*.tif'));
for i = 1:numel(tif_files)
    disp(i);
    im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));

    im = im(scale(1,1):scale(1,2),scale(2,1):scale(2,2),scale(3,1):scale(3,2));
    tifwrite(uint16(im), fullfile(write_folder, 'imgs', tif_files(i).name(1:end-4)));
end


%% crop mat files
load(fullfile(seg_res_folder_path, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d');
eig_overlay = cell(numel(eig_res_3d), 1);
for i = 1:numel(eig_res_3d)
    eig_res_2d{i} = sub_crop(eig_res_2d{i}, scale);
    eig_res_3d{i} = sub_crop(eig_res_3d{i}, scale);
%     eig_overlay{i} = sub_crop(eig_overlay{i}, scale);
%     id_mat_2nd{i} = sub_crop(id_mat_2nd{i}, scale);
end
save(fullfile(write_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d', 'eig_res_3d', 'eig_overlay', '-v7.3');
clear eig_res_2d eig_res_3d eig_overlay

load(fullfile(seg_res_folder_path, 'synQuant_refine_res_4d_v9.mat'));
for i = 1:numel(refine_res)
    refine_res{i} = sub_crop(refine_res{i}, scale);
    threshold_res{i} = sub_crop(threshold_res{i}, scale);
end
save(fullfile(write_folder, 'synQuant_refine_res_4d_v9.mat'), 'refine_res', 'threshold_res','-v7.3');
clear refine_res threshold_res

load(fullfile(seg_res_folder_path, 'varianceMap.mat'));
for i = 1:numel(varMap)
    varMap{i}{1,1} = sub_crop(varMap{i}{1,1}, scale);
    varMap{i}{1,2} = sub_crop(varMap{i}{1,2}, scale);
end
save(fullfile(write_folder, 'varianceMap.mat'), 'varMap','-v7.3');
clear varMap

load(fullfile(seg_res_folder_path, 'synQuant_res.mat'));
for i = 1:numel(fMaps)
    fMaps{i} = sub_crop(fMaps{i}, scale);
    id_mat{i} = sub_crop(id_mat{i}, scale);
    z_mat{i} = sub_crop(z_mat{i}, scale);
end
save(fullfile(write_folder, 'synQuant_res.mat'), 'fMaps', 'id_mat', 'z_mat');
clear fMaps id_mat z_mat


