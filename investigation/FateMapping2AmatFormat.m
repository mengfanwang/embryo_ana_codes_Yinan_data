% amat format:  'node id', 'type', 'x', 'y', 'z', 'radius', 'parent node id', 'timepoint', 'tag', 'lineage id'

global tracks
node_num = length(movieInfo.xCoord);
n_perframe = movieInfo.n_perframe;
tracks = zeros(node_num, 10);
tracks(:, 1) = [1:node_num]';
tracks(:, 3:5) = movieInfo.orgCoord - 1;
tracks(:, 3:4) = tracks(:,3:4) * 2;
tracks(:, 8) = movieInfo.frames - 1;

for ii = 1:node_num
    if isempty(movieInfo.parents{ii})
        tracks(ii, 7) = -1;
    else
        tracks(ii, 7) = movieInfo.parents{ii};
    end
end

% traverse: boundary cell in 1st frame
que = [];
for ii = 1:movieInfo.n_perframe(1)
    if movieInfo.zCoord(ii) >= 105
        que = [que; ii];
    end
end
while ~isempty(que)
    cur = que(1);
    tracks(cur,9) = 1;
    que = [que(2:end); movieInfo.kids{cur}];
end

% union-find: merge lineage id
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

