% clc;clear;close all;
% load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0530_000_191/movieInfo.mat');
clearvars -except movieInfo embryo_vid
load('detect_divPair.mat');

% merge devision detections to tracking results
% case disscusion:
% 1. is parent and one child in the same track? yes 1 no 2
% 2. is the child first head or the second?
movieInfo_temp = movieInfo;
for ii = 1:size(detect_divPair,1)
    if mod(ii,1000) == 0
        fprintf('%d/%d\n', ii, size(detect_divPair,1));
    end
    parent_id = detect_divPair(ii,1);
    child1_id = detect_divPair(ii,2);
    child2_id = detect_divPair(ii,3);

    parent_track = movieInfo.particle2track(parent_id,1);
    child1_track = movieInfo.particle2track(child1_id,1);
    child2_track = movieInfo.particle2track(child2_id,1);
    
    if parent_track == child1_track
        % parent connects to child1, now add to child2
        movieInfo_temp = mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child2_id, child2_track);
    else 
        % parents are not connected to any child
        movieInfo_temp = mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child1_id, child1_track);
        movieInfo_temp = mergeParent2ChildTrack...
            (movieInfo_temp, parent_id, child2_id, child2_track);
    end
end

%% save to mastodon file
mat2tgmm(movieInfo_temp, '/work/Nova/embryo_res_folder/mengfan_data_res/merge_division/tgmm_format');


function movieInfo = mergeParent2ChildTrack(movieInfo, parent_id, child_id, child_track)
    child_loc = find(movieInfo.tracks{child_track} == child_id);
    if child_loc == 1
        movieInfo.parents{child_id} = parent_id;
        movieInfo.tracks{child_track} = ...
            [parent_id; movieInfo.tracks{child_track}];
    elseif child_loc == 2
        movieInfo.parents{child_id} = parent_id;
        movieInfo.tracks{child_track} = ...
            [parent_id; movieInfo.tracks{child_track}(2:end)];
    else
        error('Wrong case');
    end
end

