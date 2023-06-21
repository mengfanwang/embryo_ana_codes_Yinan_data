clc;clear;close all;

%% get link table
link_filename = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0420/MastodonTable-Link.csv';
link_table = readtable(link_filename);

% link_table = table2array(link_table(3:end, [6 5]));
% link_ind = 1:size(link_table,1);
link_table = table2array(link_table(3:end, [4 5 7 8 9 10 11]));
link_ind = find(link_table(:,5)|link_table(:,6));
link_table = link_table(link_ind,:);

%% plot lineages
G = digraph();
EndNodes = cell(size(link_table,1),2);
for ii = 1:size(link_table,1)
    EndNodes{ii,1} = num2str(link_table(ii,2));
    EndNodes{ii,2} = num2str(link_table(ii,1));
end
G = addedge(G, EndNodes(:,1), EndNodes(:,2));
plot(G, 'Layout', 'layered');


%% get spot location
spot_filename = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0420/MastodonTable-Spot.csv';
spot_table = readtable(spot_filename);

spot_table = table2array(spot_table(3:end, [1 4 5 6 7 9 10 11 12 13])); % id time x y z
spot_ind = find(any(spot_table(:,6:10),2));
spot_table = spot_table(spot_ind,:);
% % convert to matlab corrdinate
% spot_table(:,2) = spot_table(:,2) + 1;
% spot_table(:,3:5) = spot_table(:,3:5).*[0.5 0.5 1/5.86] + 1;

% check no linked cells
inLink_flag = zeros(size(spot_table,1),1);
for ii = 1:size(spot_table)
    if ~ismember(spot_table(ii,1), link_table(:,1)) && ~ismember(spot_table(ii,1), link_table(:,2))
        inLink_flag(ii) = 1;
    end
end
% spot_table(find(inLink_flag),:) = [];

%% plot lineages in 3d
close all
spot_linked = spot_table(logical(~inLink_flag),2:5);
cmap = colormap;
for tt = 0:191
    s_ind = find(spot_linked(:,1) == tt);
    plot3(spot_linked(s_ind,2), spot_linked(s_ind,3), spot_linked(s_ind,4),'.', 'Color', cmap(193-tt,:));
    hold on
end

% get tracks
endnodes = zeros(size(EndNodes));
for ii = 1:size(endnodes, 1)
    endnodes(ii,1) = str2num(EndNodes{ii,1});
    endnodes(ii,2) = str2num(EndNodes{ii,2});
end
spot_max = max(spot_table(:,1));
spot_inverse = zeros(spot_max,1);
for ii = 1:size(spot_table,1)
    spot_inverse(spot_table(ii,1)) = ii;
end
for tt = 0:191
    s_ind = spot_table((spot_table(:,2) == tt),1);
    link_ind = endnodes((ismember(endnodes(:,1),s_ind)),:);
    s_loc = spot_table(spot_inverse(link_ind(:,1)), 3:5);
    t_loc = spot_table(spot_inverse(link_ind(:,2)), 3:5);
    x_loc = [s_loc(:,1)'; t_loc(:,1)'];
    y_loc = [s_loc(:,2)'; t_loc(:,2)'];
    z_loc = [s_loc(:,3)'; t_loc(:,3)'];
    plot3(x_loc, y_loc,z_loc, 'Color', cmap(193-tt,:));
end
    

