clc;clear;close all;
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
x = 1792;
y = 1818;
z = 253;

load('tracks_stitch_v4_1000.mat');

im = loadchannel(0);
im = im/255;
mask = loadchannel(7);  
mask = mask > 0;

% %% visual segmentations
% overlay_brain = uint8(zeros(x,y,3,z));
% for zz = 1:z
%     boundary = boundarymask(mask(:,:,zz));
%     if any(boundary, 'all')
%         overlay_brain(:,:,:,zz) = labeloverlay(im(:,:,zz),boundary,'Colormap',[1 1 0],'Transparency',0);
%     else
%         overlay_brain(:,:,:,zz) = uint8(repmat(im(:,:,zz),1,1,3));
%     end
% end

%%
tracks_last = tracks(tracks(:,8) == 1000,:);
cell_im = zeros(x,y,z);
cell_flag = zeros(size(tracks_last,1),1);

for ii = 1:size(tracks_last,1)
    cell_loc_x = round(tracks_last(ii,4));
    cell_loc_y = round(tracks_last(ii,3));
    cell_loc_z = round(tracks_last(ii,5));
    if mask(cell_loc_x, cell_loc_y, cell_loc_z)
        cell_im(cell_loc_x, cell_loc_y, cell_loc_z) = 1;
        cell_flag(ii) = 1;
    end
end

% % visualize
% cell_im = imdilate(cell_im, strel('sphere',3));
% overlay = zeros(x,y,3,z);
% for zz = 1:z
%     overlay(:,:,1,zz) = max(im(:,:,zz), cell_im(:,:,zz));
%     overlay(:,:,2,zz) = max(im(:,:,zz), cell_im(:,:,zz));
%     overlay(:,:,3,zz) = im(:,:,zz);
% end

% get linegae
select_lineage = unique(tracks_last(logical(cell_flag),10));
loc = tracks(:, 3:5);
loc(:,3) = -loc(:,3);

        writerObj = VideoWriter('myVideo.avi');
        writerObj.FrameRate = 30;
        open(writerObj);

for tt = 0:10:1000
    tt
    ind_start = find(tracks(:,8) == tt, 1, 'first');
    ind_end = find(tracks(:,8) == tt, 1, 'last');
    select_flag = false(size(tracks,1),1);
    for ii = ind_start:ind_end
        select_flag(ii) = ismember(tracks(ii,10), select_lineage);
    end

    % plot
    figure(1)
    scatter(loc(ind_start:ind_end,1), loc(ind_start:ind_end,3)*5.86+1400, 18, [0.75 0.75 0.75], 'filled');
    hold on;
    scatter(loc(select_flag,1), loc(select_flag,3)*5.86+1400, 18, [1 0 0], 'filled');
    axis([0 1800 0 1400]);
    set(gcf,'position',[100, 100, 960, 960]);

        frame = getframe(gcf);
        frame = frame.cdata;
        if tt == 0
            sz = size(frame);
            sz = sz(1:2);
        else
            frame = imresize(frame, sz);
        end
        writeVideo(writerObj, frame);
        close all;
end
        close(writerObj);



function im = loadchannel(channel)
    path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
    z = 253;
    im = zeros(1792,1818,z);
    for zz = 1:z
        im(:,:,zz) = tifread(fullfile(path, ['SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C' ...
            sprintf('%01d',channel) '_Z' sprintf('%03d',zz-1) '.tif']));
    end
end
