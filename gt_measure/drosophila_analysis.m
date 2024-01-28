a = find(fp_list);
s = RandStream('mlfg6331_64', 'Seed', 20240108);
a = randsample(s, a, 100);
% 1/2 division
kk = 1;
ii = edge(a(kk),1)
jj = edge(a(kk),2)
[movieInfo.frames(ii) movieInfo.xCoord(ii) movieInfo.yCoord(ii) movieInfo.zCoord(ii)]-1
[movieInfo.frames(jj) movieInfo.xCoord(jj) movieInfo.yCoord(jj) movieInfo.zCoord(jj)]-1

%%
im1 = tifread(fullfile(gt_detection_folder, tif_files(movieInfo.frames(ii)).name));
im2 = tifread(fullfile(gt_detection_folder, tif_files(movieInfo.frames(jj)).name));
% im1 = imdilate(im1, strel("sphere",2));
% im2 = imdilate(im2, strel("sphere",2));
data_folder = '/work/Mengfan/EmbryoData_other/drosophila-cell-tracking/raw_data/TIF';
data2 = tifread(fullfile(data_folder, [sprintf('%05d', movieInfo.frames(jj)) '.tif']));

im1(round(movieInfo.yCoord(ii)), round(movieInfo.xCoord(ii)), round(movieInfo.zCoord(ii)))
im2(round(movieInfo.yCoord(jj)), round(movieInfo.xCoord(jj)), round(movieInfo.zCoord(jj)))

%%
im2 = imdilate(im2, strel("sphere",2));
im_overlay = zeros(size(im1,1), size(im1,2), 3, size(im1,3));
for zz = 1:30
    im_overlay(:,:,:,zz) = labeloverlay(data2(:,:,zz)/500, im1(:,:,zz));
end

im1 = imdilate(im1, strel("sphere",2));
zz = round(movieInfo.zCoord(ii));
imshow(labeloverlay(data2(:,:,zz)/500, im1(:,:,zz)));

im2 = imdilate(im2, strel("sphere",2));
zz = round(movieInfo.zCoord(jj));
imshow(labeloverlay(data2(:,:,zz)/500, im2(:,:,zz)));