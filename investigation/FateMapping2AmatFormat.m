% amat format:  'node id', 'type', 'x', 'y', 'z', 'radius', 'parent node id', 'timepoint', 'tag', 'lineage id'
 
node_num = length(movieInfo.xCoord);
n_perframe = movieInfo.n_perframe;
tracks = zeros(node_num, 10);
tracks(:, 1) = [1:node_num]';
tracks(:, 3:5) = movieInfo.orgCoord - 1;
tracks(:, 8) = movieInfo.frames - 1;

%
for ii = 1:node_num
    if isempty(movieInfo.parents{ii})
        tracks(ii, 7) = -1;
    else
        tracks(ii, 7) = movieInfo.parents{ii};
    end
end

n_cumsum = [0; cumsum(n_perframe)];
for tt = 1:length(n_perframe)
    tracks(n_cumsum(tt)+1:n_cumsum(tt+1), 9) = cell_loc_all{tt}(:,4);
end

for ii = 1:length(movieInfo.tracks)
    for jj = 1:length(movieInfo.tracks{ii})
        tracks(movieInfo.tracks{ii}(jj),10) = ii;
    end
end