clc
% plot a few tracks, especially in the split-merge module

%% given time and a voxel, find corresponing cell and track
tt = 12;
data_size = size(refine_res{tt});

% loc = [70 41 33];
loc = [67 44 33];
% loc = [66 52 33];
query = sub2ind(data_size, loc(1), loc(2), loc(3));
for ii = 1:length(movieInfo.xCoord)
    if movieInfo.frames(ii) == tt
        if ismember(query, movieInfo.voxIdx{ii})
            point = ii;
            break
        end
    end
end
for ii = 1:length(movieInfo.tracks)
    if ismember(point, movieInfo.tracks{ii})
        track = ii
    end
end
fprintf('Point: %d Track: %d\n', point, track);

%% find all tracks before merge
a = movieInfo_noJump.tracks{track};
b = cellfun(@(x)any(ismember(a,x)), movieInfo_noJump.track_bf_merge, 'UniformOutput', false);
c = [];
for ii = 1:length(b);if b{ii}; c = [c ii];end;end
tracks = movieInfo_noJump.track_bf_merge(c);

%% plot all given tracks
G = digraph();
for ii = 1:length(tracks)
    for jj = 1:length(tracks{ii})-1
        G = addedge(G,num2str(tracks{ii}(jj)), num2str(tracks{ii}(jj+1)));
    end
end
plot(G, 'Layout', 'layered');

