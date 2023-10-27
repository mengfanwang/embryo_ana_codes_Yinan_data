clc
% given the cell location, find the cell id and the track id

% movieInfo = movieInfoAll{1};
query = [1414.8 562.2 1];
[ii, track] = findPoint(query, movieInfo);


% query = [1038.1 317.8 73.4];


% query = ([396.3960  325.4381  145.6634]-1).*[2 2 1];
% [jj, track] = findPoint(query, movieInfo);

function [point, track] = findPoint(query, movieInfo)
    query = query.*[0.5 0.5 1] + 1;
%     query = query.*[0.5 0.5 1/5.86] + 1;
    first_min = inf;
    second_min = inf;
    for ii = 1:length(movieInfo.xCoord)
%         current_dist = norm(query - movieInfo.orgCoord(ii,:));
        current_dist = norm(query - [movieInfo.xCoord(ii) movieInfo.yCoord(ii) movieInfo.zCoord(ii)]);
        if current_dist < first_min
            second_min = first_min;
            first_min = current_dist;
            point = ii;
        elseif current_dist < second_min
            second_min = current_dist;
        end 
    end
    for ii = 1:length(movieInfo.tracks)
        if ismember(point, movieInfo.tracks{ii})
            track = ii;
        end
    end
    fprintf('Point: %d Track: %d\n', point, track);
    fprintf('Distance: %f 2nd distance: %f\n', first_min, second_min);
end