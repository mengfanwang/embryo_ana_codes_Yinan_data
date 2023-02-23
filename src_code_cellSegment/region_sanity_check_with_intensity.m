function regLabelOut = region_sanity_check_with_intensity(vid, regLabel, minSz, minIntensity)


ids = unique(regLabel(regLabel>0));

stats_new = regionprops3(regLabel,'Volume', 'VoxelIdxList');
regCnt = 0;
   regLabelOut = zeros(size(regLabel));
for i=1:length(ids)
    voxIdx = stats_new.VoxelIdxList{ids(i)};
    cur_mean_intensity = mean(double(vid(voxIdx)));
    if length(voxIdx) <= minSz || cur_mean_intensity <= minIntensity
        continue;
    end
    regCnt = regCnt + 1;
    regLabelOut(voxIdx) = regCnt;
end