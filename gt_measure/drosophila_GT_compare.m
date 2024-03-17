clc;clear;close all;

file_path = '/work/public/Embryo_annotation';
resolution = [1 1 1];
radius = 30;

order = {'Public', 'Zeyuan', 'Kevin', 'Our detection'};
node = cell(4,1);
edge = cell(4,1);

%% public result
node_file = fullfile(file_path, 'MastodonTable_public-Spot.csv');
edge_file = fullfile(file_path, 'MastodonTable_public-Link.csv');
% [node{1}, edge{1}] = readMastodon(node_file, edge_file);
[node_public, edge_public] = readMastodon(node_file, edge_file);

%% zeyuan result
node_file = fullfile(file_path, 'MastodonTable_zeyuan-Spot.csv');
edge_file = fullfile(file_path, 'MastodonTable_zeyuan-Link.csv');
% [node{2}, edge{2}] = readMastodon(node_file, edge_file);
[node_zeyuan, edge_zeyuan] = readMastodon(node_file, edge_file);

%% kevin result
node_file = fullfile(file_path, 'MastodonTable_kevin-Spot.csv');
edge_file = fullfile(file_path, 'MastodonTable_kevin-Link.csv');
% [node{3}, edge{3}] = readMastodon(node_file, edge_file);
[node_kevin, edge_kevin] = readMastodon(node_file, edge_file);

%% load our detection result
path = '/work/Mengfan/EmbryoData_other/drosophila-cell-tracking/Our/Tracking/';
files = dir(path); files = files(3:end-1);
% load('/work/Mengfan/EmbryoData_other/drosophila-cell-tracking/Our/Tracking/0113_iter10_0_99_z05s3/movieInfo.mat');
% node:[id x y z t]
diary('log.txt');
for ff = 10:10
    load(fullfile(path, files(ff).name, 'movieInfo.mat'));
node_num = length(movieInfo.xCoord);
node = zeros(node_num, 5);
node(:,1) = (1:node_num)';
node(:,2:4) = movieInfo.orgCoord;
node(:,5) = movieInfo.frames;

edge = zeros(node_num, 2);
edge_num = 0;
for ii = 1:length(movieInfo.tracks)
    for jj = 1:length(movieInfo.tracks{ii})-1
        edge_num = edge_num + 1;
        edge(edge_num,:) = [movieInfo.tracks{ii}(jj) movieInfo.tracks{ii}(jj+1)]; 
    end
end
edge = edge(1:edge_num, :);

[tp, fp, fn] = evaluation(node_public, edge_public, ...
            [node radius*ones(size(node,1),1)], edge, resolution);
        fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f\n', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));
end
%% compare
for ii = 1:3
    for jj = 1:3
        if ii ~= jj
        [tp, fp, fn] = evaluation(node{ii}, edge{ii}, ...
            [node{jj} radius*ones(size(node{jj},1),1)], edge{jj}, resolution);
        fprintf('%s to %s:\n', order{ii}, order{jj});
        fprintf('Precision: %2.2f Recall: %2.2f Accuracy: %2.2f\n', tp*100/(tp+fp), tp*100/(tp+fn), tp*100/(tp+fp+fn));
        end
    end
end


function [node, edge] = readMastodon(node_file, edge_file)
    % node:[id x y z t]
    % edge:[id 1 id2]
    node = readmatrix(node_file);
    node = node(3:end, [2 5 6 7 4]);
    node(:,1) = node(:,1) + 1;
    edge = readmatrix(edge_file);
    edge = edge(3:end, [4 5]);
    if edge(1,1) > edge(1,2)
        edge = edge(:,[2 1]);
    end
    edge = edge + 1;
end