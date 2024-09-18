    clc;clear;close all;
% plot the division heatmap
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
x = 1792;
y = 1818;
z = 253;

load('../tracks_stitch_v4_1000.mat');
global tracks
cell_num = size(tracks,1);
loc = tracks(:, 3:5);
loc(:,3) = -loc(:,3);
tracks(:, 10) = [1:cell_num]';

    
channel_list = 0:8;
name = {'All', 'Right Eye', 'Left Eye', 'Brain', 'Somite AB', 'Somite BA', 'Somite AA', 'Somite BB', 'Tail bud'};
color = 1-(1-[0, 447/1000, 741/1000]).*[0.1:0.9/127:1]';

% get division num
child_cnt = zeros(cell_num,1);
for ii = 1:cell_num
    if tracks(ii,7) > 0
        child_cnt(tracks(ii,7)) = child_cnt(tracks(ii,7)) + 1;
    end
end



cnt = 101;
for tt = 1000:-10:0
    tt
    ind_start = find(tracks(:,8) == tt, 1, 'first');
    ind_end = find(tracks(:,8) == tt, 1, 'last');

    for ii = ind_start:cell_num
        if tracks(ii,7) > 0
            addChild(tracks(ii,7), ii);
        end
    end
    for ii = ind_start:cell_num
        tracks(ii,10) = findRoot(tracks(ii,10));
    end

    division_cnt = zeros(cell_num,1);
    for ii = ind_start:cell_num
        if child_cnt(ii) == 2
            division_cnt(tracks(ii,10)) = division_cnt(tracks(ii,10)) + 1;
        elseif child_cnt(ii) > 2
            error('More than two kids!');
        end
    end

    
%     select_flag = false(size(tracks,1),1);
%     for ii = ind_start:ind_end
%         select_flag(ii) = ismember(tracks(ii,10), select_lineage);
%     end
    select_flag = unique(tracks(ind_start:ind_end,10));

    % plot
    colormap(color);
    scatter(loc(select_flag,1), loc(select_flag,3)*5.86+1400, 18, min(division_cnt(select_flag),...
        quantile(division_cnt(select_flag), 0.98)), 'filled');
    axis([0 1800 0 1400]);
    set(gcf,'position',[100, 100, 960, 960]);
    text(50, 50, sprintf('T = %d', tt), 'FontSize', 18);
    c = colorbar;

legend(name{channel_list+1});
        frame = getframe(gcf);
        frame = frame.cdata;
        if tt == 1000
            sz = size(frame);
            sz = sz(1:2);
            video = zeros(sz(1), sz(2), 3, 101);
        else
            frame = imresize(frame, sz);
        end
        video(:,:,:,cnt) = frame;
        cnt = cnt - 1;

end


        writerObj = VideoWriter('Division.avi');
        writerObj.FrameRate = 30;
        open(writerObj);

        for tt = 1:101
            writeVideo(writerObj, uint8(video(:,:,:,tt)));
        end
        close(writerObj);

function root = findRoot(u)
    global tracks
    if tracks(u, 7) > 0 && tracks(u,10) ~= u
        root = findRoot(tracks(u,10));
        tracks(u,10) = root;
    else
        root = u;
    end
end

function addChild(u, v)
    % u is the parent of v
    global tracks
    root_u = findRoot(u);
    root_v = findRoot(v);
    tracks(root_v, 10) = root_u;
end

function im = loadchannel(channel)
    path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
    z = 253;
    im = zeros(1792,1818,z);
    for zz = 1:z
        im(:,:,zz) = tifread(fullfile(path, ['SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C' ...
            sprintf('%01d',channel) '_Z' sprintf('%03d',zz-1) '.tif']));
    end
end
