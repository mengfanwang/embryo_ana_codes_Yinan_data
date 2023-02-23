% crop a small portion of the embryo data
save_folder = '/home/congchao/Desktop/wcc_embryo/yinan_data/crop_data_scale_0.5_temporal_aligned';

scale = [200 450; 200 450; 125, 175];
data_folder = '/home/congchao/Desktop/wcc_embryo/yinan_data/full_data_scale_0.5_temporal_aligned';

tif_files = dir(fullfile(data_folder, '*.tif'));

for i=1:numel(tif_files)
    disp(i);
    im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    im = im(scale(1,1):scale(1,2),scale(2,1):scale(2,2),scale(3,1):scale(3,2));
    tifwrite(uint16(im), ...
        fullfile(save_folder, tif_files(i).name));
end

