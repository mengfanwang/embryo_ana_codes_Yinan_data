% movieInfo = m191.movieInfo;
movieInfo = movieInfoAll{15};
                                

%% check jump num
jumpNum = zeros(5,1);
for ii = 1:length(movieInfo.tracks)
    jumpTmp = diff(movieInfo.frames(movieInfo.tracks{ii}));
    for jj = 1:length(jumpTmp)
        jumpNum(jumpTmp(jj)) = jumpNum(jumpTmp(jj)) + 1;
    end
end
% jumpNum = jumpNum/sum(jumpNum);

% jumpNum

1-jumpNum(1)

%% check avg track length without allowing jump
track_length = cell(length(movieInfo.tracks), 1);
for ii = 1:length(movieInfo.tracks)
    jumpTmp = diff(movieInfo.frames(movieInfo.tracks{ii}));
    jumpPt = find(jumpTmp>1);       
    if isempty(jumpPt)
        track_length{ii} = length(movieInfo.tracks{ii});
    elseif length(jumpPt) == 1
        track_length{ii} = [jumpPt length(movieInfo.tracks{ii})-jumpPt];
    else
        track_length{ii} = jumpPt(1);
        for jj = 1:length(jumpPt)-1
            track_length{ii} = [track_length{ii} jumpPt(jj+1)-jumpPt(jj)];
        end
        track_length{ii} = [track_length{ii} length(movieInfo.tracks{ii})-jumpPt(end)];
    end
end
track_num = 0;
track_sum = 0;
for ii = 1:length(movieInfo.tracks) 
    track_num = track_num + length(track_length{ii});
    track_sum = track_sum + sum(track_length{ii});
end
track_sum/track_num
