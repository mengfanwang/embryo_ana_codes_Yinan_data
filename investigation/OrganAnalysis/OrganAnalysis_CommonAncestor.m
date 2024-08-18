clc;clear;close all;
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
% SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C0_Z000.tif
x = 1792;
y = 1818;
z = 253;

load('../tracks_stitch_v4_1000.mat');

channel_list = 1:8;
name = {'Right Eye', 'Left Eye', 'Brain', 'Somite AB', 'Somite BA', 'Somite AA', 'Somite BB', 'Tail bud'};

tracks_last = tracks(tracks(:,8) == 1000,:);
cell_flag_list = zeros(size(tracks_last,1),length(channel_list));
cell_im = zeros(x,y,z);
for cc = 1:length(channel_list)
    mask = loadchannel(channel_list(cc));  
    mask = mask > 0;
   
    if channel_list(cc) == 0
        cell_flag = ones(size(tracks_last,1),1);
    else
        cell_flag = zeros(size(tracks_last,1),1);    
        for ii = 1:size(tracks_last,1)
            cell_loc_x = round(tracks_last(ii,4));
            cell_loc_y = round(tracks_last(ii,3));
            cell_loc_z = round(tracks_last(ii,5));
            if mask(cell_loc_x, cell_loc_y, cell_loc_z)
                if cell_im(cell_loc_x, cell_loc_y, cell_loc_z) > 0
                    warning('Cell has mutilple labels!');
                end
                cell_im(cell_loc_x, cell_loc_y, cell_loc_z) = cc;
                cell_flag(ii) = 1;
            end
        end
    end
    cell_flag_list(:,cc) = cell_flag;
end


for tt = 0:0
    ind_start = find(tracks(:,8) == tt, 1, 'first');
    ind_end = find(tracks(:,8) == tt, 1, 'last');
    select_flag = zeros(ind_end-ind_start+1,1);                             % 8189*1
    select_tag = zeros(ind_end-ind_start+1,length(channel_list));           % 8189*8
for cc = 1:length(channel_list)
    % get linegae
    select_lineage = unique(tracks_last(logical(cell_flag_list(:,cc)),10));
    
    for ii = ind_start:ind_end
        if ismember(tracks(ii,10), select_lineage)
            select_tag(ii,cc) = 1;
        
            if select_flag(ii) > 0
                warning('Common ancestor!');
            else
                select_flag(ii) = cc;
            end
        end
    end
end
end
common_ancestor = find(sum(select_tag,2)>1);
varNames = {'Node', 'Initial Location', 'Tag', 'End Location'};
end_location = cell(size(common_ancestor,1),1);
for ii = 1:size(common_ancestor,1)
    end_location{ii} = tracks_last(tracks_last(:,10) == tracks(common_ancestor(ii),10) & any(cell_flag_list,2), 3:5);
end
t = table(common_ancestor, tracks(common_ancestor,3:5), ...
        array2table(select_tag(common_ancestor,:), 'VariableNames', name), ...
        end_location, 'VariableNames', varNames);
% find error at tracks_last(5118). its a somite AB but not shown in
% mastodon. location: [1171.43362831858,835.840707964602,217.982300884956]

function im = loadchannel(channel)
    path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
    z = 253;
    im = zeros(1792,1818,z);
    for zz = 1:z
        im(:,:,zz) = tifread(fullfile(path, ['SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C' ...
            sprintf('%01d',channel) '_Z' sprintf('%03d',zz-1) '.tif']));
    end
end
