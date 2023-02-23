% clc
% movieInfo = movieInfoAll{23};


query = [1204.6 704.2 92.1];
[ii, track] = findPoint(query, movieInfo);

query = [1181.3 702.5 95];
[jj, ~] = findPoint(query, movieInfo);
% movieInfo.tracks{track} 
% movieInfo.nei{point}
% movieInfo.frames(movieInfo.nei{point})
% movieInfo.CDist{point}

function [point, track] = findPoint(query, movieInfo)
    query = query.*[0.5 0.5 1] + 1;
    first_min = inf;
    second_min = inf;
    for ii = 1:length(movieInfo.xCoord)
        current_dist = norm(query - movieInfo.orgCoord(ii,:));
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