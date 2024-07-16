clc;clear;
% convert the tgmm result of Amat2014 data to n*10 column format
% 'node id', 'type', 'x', 'y', 'z', 'radius', 'parent node id', 'timepoint', 'tag', 'lineage id'

%% load data
% file_path = '/work/Mengfan/TGMM2.0/Outputs/Tracking_results/GMEMtracking3D_2024_3_6_20_13_21/';
% node_file = fullfile(file_path, 'MastodonTable-Spot.csv');
% edge_file = fullfile(file_path, 'MastodonTable-Link.csv');
node_file = '/work/Nova/embryo_res_folder/mengfan_data_res/ground_truth/Elephant_detection_410_Spot.csv';
edge_file = '/work/Nova/embryo_res_folder/mengfan_data_res/ground_truth/MastodonTable-0_192-Link-new-Link.csv';

node = readmatrix(node_file);
node = node(3:end, [2 5 6 7 4]);
edge = readmatrix(edge_file);
edge = edge(3:end, [4 5]);
edge = edge + 1;

%% ultrack
ul_file = '/work/Mengfan/EmbryoData_other/ultrack_tracking_results/Yinan.csv';
ul = readmatrix(ul_file);
node = ul(:,5:-1:2);
node = [ul(:,6),node(:,1:end)];
[~, ind] = sort(node(:,5));
node = node(ind,:);

edge = [ul(:,8),ul(:,6)];
ul_id = find(edge(:,1)==-1);
edge(ul_id,:)=[];

for ii = 1:size(edge,1)
    for jj = 1:2
        edge(ii,jj) = find(node(:,1) == edge(ii,jj));
    end
end

%% convert to n*10 format
global tracks
node_num = size(node,1);
tracks = zeros(node_num, 10);
tracks(:, 1) = [1:node_num]';
tracks(:, 3:5) = node(:,2:4);
tracks(:,8) = node(:,5);

% get parent
tracks(:,7) = -1;
for ii = 1:size(edge,1)
    tracks(edge(ii,2),7) = edge(ii,1);
end

% union-find: merge lineage id
fprintf('Merge lineage id...');

tracks(:, 10) = [1:node_num]';
tic;
for ii = 1:node_num
    if tracks(ii,7) ~= -1
        addChild(tracks(ii,7), ii);
    end
end
for ii = 1:node_num
    tracks(ii,10) = findRoot(tracks(ii,10));
end
toc


% use tag == 2 to label new tracks
for ii = 1:node_num
    if tracks(tracks(ii,10), 8) > 1
        tracks(ii,9) = 2;
    end
end

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

