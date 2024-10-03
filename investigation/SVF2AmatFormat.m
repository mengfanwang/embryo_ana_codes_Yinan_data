clc;clear;

file_name = '/run/user/1008/gvfs/smb-share:server=rs0001.local,share=shared_yinan/Amat2014/SVF/SVFDatabase-SVF.csv';
file = readmatrix(file_name);

global tracks
node_num = size(file,1);
tracks = zeros(node_num, 10);
tracks(:, 1) = [1:node_num]';
tracks(:, 3:5) = file(:,3:5);
tracks(:, 7) = -1;
tracks(:, 8) = file(:,9);

for ii = 1:node_num
    if file(ii, 2) ~= -1
        tracks(ii, 7) = find(file(ii, 2) == file(:,1), 1);
    end
end



fprintf('Merge lineage id...');

tracks(:, 10) = [1:node_num]';
tic;
for ii = 1:node_num
    if ~isempty(movieInfo.parents{ii})
        addChild(movieInfo.parents{ii}, ii);
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

