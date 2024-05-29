clearvars -except movieInfo

close all;
c = colormap('turbo');
t = length(movieInfo.n_perframe);
figure(1);hold on;
cnt = 0;
matching = cell(1,t);
for ii = 1:1000 %length(movieInfo.tracks)
    if mod(ii, 1000) == 0
        fprintf("%d/%d\n", ii, length(movieInfo.tracks));
    end
    if length(movieInfo.tracks{ii}) < 10
        continue
    end
    loc = movieInfo.orgCoord(movieInfo.tracks{ii}([1 end]),:);
    if norm(loc(1,:)-loc(2,:)) < 100
        continue
    end
    for jj = 1:length(movieInfo.tracks{ii})-1
        tt = movieInfo.frames(movieInfo.tracks{ii}(jj));
        loc = movieInfo.orgCoord(movieInfo.tracks{ii}(jj:jj+1),:);
        matching{tt} = [matching{tt}; loc(1,:); loc(2,:); nan(1,3);];
    end
    cnt = cnt + 1;
end

%%
for tt = 1:191
    tt
    loc = matching{tt};
    plot(loc(:,1), loc(:,2), 'Color', c(round(tt/t*256),:));
%     plot3(loc(:,1), loc(:,2), loc(:,3)*5.4, 'Color', c(round(tt/t*256),:));
end
axis([-100 900 0 1000]);
set(gca, 'Color', 'none'); 