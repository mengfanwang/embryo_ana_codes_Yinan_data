clc;clear;close all;

%% get detection result
im_size = [1920 1920 180];
z_threshod = 0.5;
voxIdx = cell(2e6,1);
vox_num = 0;
file_name = '/work/Mengfan/TGMM2.0/Outputs/Tracking_results/21-01-11/XML_finalResult_lht/GMEMfinalResult_frame';
n_perframe = zeros(192,1);
for tt = 0:191
    tt
    xml_data = xml2struct([file_name num2str(tt, '%.4d') '.xml']);
    xml_data = xml_data.document.GaussianMixtureModel;
    idmap = zeros(im_size);
    for ii = 1:length(xml_data)
        mu = str2num(xml_data{ii}.Attributes.m);
        preMatrix = str2num(xml_data{ii}.Attributes.W);
        preMatrix = reshape(preMatrix, [3 3]);
        covMatrix = inv(preMatrix);
        
        sigma = sqrt(sum(abs(covMatrix)));
        low_bound = round(max(mu - 2*sigma, 0));
        up_bound = round(min(mu + 2*sigma, im_size-1));
        % xyz corrdinate are real world corrdniate
        [x_corrd, y_corrd, z_corrd] = meshgrid(low_bound(1):up_bound(1), low_bound(2):up_bound(2), low_bound(3):up_bound(3));
        xyz_corrd = cat(4, x_corrd, y_corrd, z_corrd);
        xyz_corrd = reshape(xyz_corrd, [numel(x_corrd) 3]);
        mahal_dist = sqrt(sum((xyz_corrd-mu)*preMatrix.*(xyz_corrd-mu), 2));
        

        valid_corrd = find(mahal_dist < z_threshod);
        idmap_temp = idmap(low_bound(2)+1:up_bound(2)+1, low_bound(1)+1:up_bound(1)+1, low_bound(3)+1:up_bound(3)+1);
        if isempty(valid_corrd)  % at least give one pixel
            valid_corrd = sub2ind(size(idmap_temp), round(mu(2))-low_bound(2)+1, round(mu(1))-low_bound(1)+1, round(mu(3))-low_bound(3)+1);
        end
        idmap_temp(valid_corrd) = ii;
        idmap(low_bound(2)+1:up_bound(2)+1, low_bound(1)+1:up_bound(1)+1, low_bound(3)+1:up_bound(3)+1) = idmap_temp;
        
        vox_num = vox_num + 1;
    end
    n_perframe(tt+1) = vox_num;
    voxIdx_temp = regionprops3(idmap, 'VoxelIdxList');
    voxIdx_temp = voxIdx_temp.VoxelIdxList;
    if tt == 0
        voxIdx(1:vox_num) = voxIdx_temp;
    else
        voxIdx(n_perframe(tt)+1:vox_num) = voxIdx_temp;
    end    
end
voxIdx = voxIdx(1:vox_num);


