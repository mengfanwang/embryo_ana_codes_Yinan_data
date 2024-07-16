clc;clear;close all;
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
z = 253;

load('tracks_stitch_v4_1000.mat');

im = loadchannel(0);
im = im/255;
mask_brain = loadchannel(3);
mask_brain = mask_brain > 0;

% %% visual segmentations
% overlay_brain = uint8(zeros(1792,1818,3,z));
% for zz = 1:z
%     boundary = boundarymask(mask_brain(:,:,zz));
%     if any(boundary, 'all')
%         overlay_brain(:,:,:,zz) = labeloverlay(im(:,:,zz),boundary,'Colormap',[1 1 0],'Transparency',0);
%     else
%         overlay_brain(:,:,:,zz) = uint8(repmat(im(:,:,zz),1,1,3));
%     end
% end

%%
tracks_last = tracks(tracks(:,8) == 1000,:);




function im = loadchannel(channel)
    path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
    z = 253;
    im = zeros(1792,1818,z);
    for zz = 1:z
        im(:,:,zz) = tifread(fullfile(path, ['SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C' ...
            sprintf('%01d',channel) '_Z' sprintf('%03d',zz-1) '.tif']));
    end
end
