clc;clear;close all;
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
x = 1792;
y = 1818;
z = 253;

load('tracks_stitch_v4_1000.mat');

channel_list = 0:8;
name = {'All', 'Right Eye', 'Left Eye', 'Brain', 'Somite AB', 'Somite BA', 'Somite AA', 'Somite BB', 'Tail bud'};

% color = [0.75 0.75 0.75; 
%          0 0.4470 0.7410;
%          0.8500 0.3250 0.0980;
%          0.9290 0.6940 0.1250;
%          0.4940 0.1840 0.5560;
%          0.4660 0.6740 0.1880;
%          0.3010 0.7450 0.9330;
%          0.6350 0.0780 0.1840];

color = [0.75 0.75 0.75;
         1    0    0;
         0    1    0;
         0    0    1;
         0    1    1;
         1    0    1;
         0.9290 0.6940 0.1250;
         0.4940 0.1840 0.5560;
         0.6350 0.0780 0.1840];

cell_flag_list = cell(length(channel_list),1);
for cc = 1:length(channel_list)
    mask = loadchannel(channel_list(cc));  
    mask = mask > 0;
    
    tracks_last = tracks(tracks(:,8) == 1000,:);
    cell_im = zeros(x,y,z);

    if channel_list(cc) == 0
        cell_flag = ones(size(tracks_last,1),1);
    else
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
    end
    cell_flag_list{cc} = cell_flag;
end

        writerObj = VideoWriter('myVideo.avi');
        writerObj.FrameRate = 30;
        open(writerObj);

for tt = 0:10:1000

for cc = 1:length(channel_list)
    % get linegae
    select_lineage = unique(tracks_last(logical(cell_flag_list{cc}),10));
    loc = tracks(:, 3:5);
    loc(:,3) = -loc(:,3);

    ind_start = find(tracks(:,8) == tt, 1, 'first');
    ind_end = find(tracks(:,8) == tt, 1, 'last');
    select_flag = false(size(tracks,1),1);
    for ii = ind_start:ind_end
        select_flag(ii) = ismember(tracks(ii,10), select_lineage);
    end

    % plot
    scatter(loc(select_flag,1), loc(select_flag,3)*5.86+1400, 18, color(cc,:), 'filled');
    hold on;
%     axis([0 1800 400 1400]);
%     set(gcf,'position',[100, 100, 960, 700]);
    axis([0 1800 0 1400]);
    set(gcf,'position',[100, 100, 960, 960]);
end
legend(name{channel_list+1});
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
