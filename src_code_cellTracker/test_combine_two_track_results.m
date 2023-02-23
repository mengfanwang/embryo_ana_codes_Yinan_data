% test merge two batches of tracking results
% only works on cbil server.
base_dir = '/home/ccwang/Desktop/embryo_res_folder/mengfan_data_res/';
batch1_folder = 'mengfan_single_view_50_69_10072022';
batch2_folder = 'mengfan_single_view_66_88_10072022';
m1 = load(fullfile(base_dir, batch1_folder, 'movieInfo.mat'));
r1 = load(fullfile(base_dir, batch1_folder, 'refine_res.mat'));
m2 = load(fullfile(base_dir, batch2_folder, 'movieInfo.mat'));

ov = 4;
movieInfo = combine_movieInfo(m1.movieInfo, m2.movieInfo, r1.out_refine_res, ov);

save_folder = '/home/ccwang/Desktop/embryo_res_folder/mengfan_data_res/combined_50-88_10162022';
data_folder = '/work/public/sameViewFusion/merged_batch/sameViewFusion_050-088_11';
mastodon_dir = fullfile(save_folder, 'mastodon');
if ~exist(mastodon_dir)
    mkdir(mastodon_dir);
end
addpath('TGMM_wrapper/');
mat2tgmm(movieInfo, fullfile(mastodon_dir, 'tgmm_format'));
tif2bdv(data_folder, fullfile(mastodon_dir, 'embryo_data_h5'));