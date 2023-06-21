clc
% given the cell location, find the cell id and the track id

% movieInfo = movieInfoAll{1};


query = [119.8 288.3 17];
[ii, track] = findPoint(query, movieInfo);

query = [117.5 271.0 17.3];
[jj, track] = findPoint(query, movieInfo);

% query = [193.78260000000	231.260900000000	16.500000000000];
% [ii, track] = findPoint(query, movieInfo);
% query = [175.9 78.2 3.1];
% query = [1099.92020000000	621.818800000000	723.073018000000];
% [jj, track2] = findPoint(query, movieInfo);
% movieInfo.tracks{track} 
% movieInfo.nei{point}[
% movieInfo.frames(movieInfo.nei{point})
% movieInfo.CDist{point}

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