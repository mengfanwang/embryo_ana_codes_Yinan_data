clc;clear;close all;
% plot the division heatmap
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
x = 1792;
y = 1818;
z = 253;


channel_list = [1 2 3];
name = {'All', 'Right Eye', 'Left Eye', 'Brain', 'Somite AB', 'Somite BA', 'Somite AA', 'Somite BB', 'Tail bud'};


load('../tracks_stitch_v4_1000.mat');
global tracks
cell_num = size(tracks,1);
loc = tracks(:, 3:5);
loc(:,3) = -loc(:,3);

load CommonAncestor.mat
load tracks_back.mat
common_list = [];
for ii = 1:size(t,1)
    if any(table2array(t.Tag(ii,channel_list)))
        common_list = [common_list; find(tracks(:,10) == tracks(t.Node(ii),10))];
    end
end

tracks(:, 10) = [1:cell_num]';

    


color = [1 0 0;
         0 1 0;
         0 0 1;];
color = 1 - 2*(1 - color)/3;

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


%%
tic;
cnt = 101;
% tracks_back = cell(101,1);
for tt = 1000:-10:0
    tt

    ind_start = find(tracks(:,8) == tt, 1, 'first');
    ind_end = find(tracks(:,8) == tt, 1, 'last');

%     for ii = ind_start:cell_num
%         if tracks(ii,7) > 0
%             addChild(tracks(ii,7), ii);
%         end
%     end
%     for ii = ind_start:cell_num
%         tracks(ii,10) = findRoot(tracks(ii,10));
%     end
%     tracks_back{cnt} = tracks(:,10);
    tracks(:,10) = tracks_back{cnt};
    tracks_last = tracks(tracks(:,8) == 1000,:);

common_flag = false(size(tracks,1),1);
for cc = 1:length(channel_list)
    select_lineage = unique(tracks_last(logical(cell_flag_list{cc}),10));
    select_flag = false(size(tracks,1),1);
    for ii = ind_start:ind_end
        select_flag(ii) = ismember(tracks(ii,10), select_lineage);
    end
    common_flag = common_flag | select_flag;

    % plot
    scatter(loc(select_flag,1), loc(select_flag,3)*5.86+1400, 18, color(cc,:), 'filled');
    hold on;
end
for ii = ind_start:ind_end
    if common_flag(ii)
        common_flag(ii) = ismember(ii, common_list);
    end
end

scatter(loc(common_flag,1), loc(common_flag,3)*5.86+1400, 30, [0 0 0], 'filled');
axis([0 1800 0 1400]);
set(gcf,'position',[100, 100, 960, 960]);
text(50, 50, sprintf('T = %d', tt), 'FontSize', 18);
legend(name{channel_list+1},'Common ancestor');
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
        close all;

end
toc

        writerObj = VideoWriter('CommonAncestor.avi');
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
