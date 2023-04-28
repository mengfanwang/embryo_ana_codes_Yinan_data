clc;

% given two different results, find all
% different tracks

% cells distance < 1 is regarded as the same cell

%%
movieInfo1 = movieInfo_old;
movieInfo2 = movieInfo_new;
paired_track_in1 = findPairedTrack(movieInfo1, movieInfo2);
paired_track_in2 = findPairedTrack(movieInfo2, movieInfo1);

%% write results
movieInfo1.tracks = movieInfo1.tracks(paired_track_in2<0);
mat2tgmm(movieInfo1, fullfile(mastodon_dir, 'ComparedResult_old'));

movieInfo2.tracks = movieInfo2.tracks(paired_track_in1<0);
mat2tgmm(movieInfo2, fullfile(mastodon_dir, 'ComparedResult_new'));


%%
function paired_track = findPairedTrack(movieInfo1, movieInfo2) 
    tracks1 = movieInfo1.tracks;
    tracks2 = movieInfo2.tracks;
    % find tracks2 in tracks1
    paired_track = -1*ones(length(tracks2),1);
    for ii = 1:length(tracks2)
        track2 = tracks2{ii};
        flag = false;
        for jj = 1:length(tracks1)
            track1 = tracks1{jj};
            % length different
            if length(track1) ~= length(track2)
                continue;
            end
            if ~all(movieInfo1.frames(track1) == movieInfo2.frames(track2))
                continue;
            end
            loc1 = movieInfo1.orgCoord(track1,:)';
            loc2 = movieInfo2.orgCoord(track2,:)';
            track_dist = vecnorm(loc1-loc2);
            if all(track_dist < 1)
                flag = true;
                paired_track(ii) = jj;
                break;
            end
        end
    end
end