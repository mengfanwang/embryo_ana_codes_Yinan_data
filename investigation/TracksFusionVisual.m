%% load data
% files = {'view9_0218_0_19', 'view10_0218_0_19', 'view11_0221_0_19', 'view12_0214_0_19'};
files = {'view9_0_30', 'view10_0_30', 'view11_0_30', 'view12_0_30'};
path = '/work/Mengfan/Embryo/20240324 isl2bGFP H2A mCherry 30tp test/Tracking';
movieInfoAll = cell(4,1);
for ff = 1:4
    load(fullfile(path, files{ff}, 'movieInfo.mat'));
    movieInfoAll{ff} = movieInfo;
end

%% load transformation matrix
path = '/work/Mengfan/Embryo/20240324 isl2bGFP H2A mCherry 30tp test';
xml_name = fullfile(path, '20240312 isl2bGFP H2AmCherry.xml');
addpath ./gt_measure/ 

t = 20;
num_view = 16;
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;

trans = cell(t, 1);
for tt = 1:t
    trans{tt} = cell(4,1);
    for vv = 1:4
        ind = (tt-1)*num_view + (vv+8);
        trans{tt}{vv} = eye(4,4);
        for jj = 1:length(register_info{ind}.ViewTransform)
            trans_temp = cellfun(@str2num, split(register_info{ind}.ViewTransform{jj}.affine.Text));
            trans_temp = [reshape(trans_temp,4,3)'; 0 0 0 1];
            trans{tt}{vv} = trans{tt}{vv}*trans_temp;
        end
%         trans{vv}(1:3,:) = trans{vv}(1:3,:) / downsample_scale;  % downsample
    end
end
%%
clearvars -except movieInfoAll trans
close all;
c = colormap('turbo');
matchingAll = cell(1,4);
for ff = 1:4
    movieInfo = movieInfoAll{ff};
    t = length(movieInfo.n_perframe);
    figure(ff);hold on;
    cnt = 0;
    matching = cell(1,t);
    for ii = 1:length(movieInfo.tracks)
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
            loc = loc(:, [2 1 3]);
            loc = (loc - 1) .* [2 2 1];
            loc1 = trans{tt}{ff} * [loc(1,:) 1]';
            loc1 = loc1(1:3)';
            loc2 = trans{tt+1}{ff} * [loc(2,:) 1]';
            loc2 = loc2(1:3)';
            matching{tt} = [matching{tt}; loc1; loc2; nan(1,3);];
        end
        cnt = cnt + 1;
    end
    matchingAll{ff} = matching;
    
    for tt = 1:t-1
        tt
        loc = matching{tt};
        plot(loc(:,2), loc(:,3), 'Color', c(round(tt/t*256),:));
    %     plot3(loc(:,1), loc(:,2), loc(:,3)*5.4, 'Color', c(round(tt/t*256),:));
    end
    % axis([-100 900 0 1000]);
    set(gca, 'Color', 'none'); 
end

%% overallfig
figure(5); hold on;
for tt = 1:t-1
    tt
    for ff = 1:4
        loc = matchingAll{ff}{tt};
        plot(loc(:,2), loc(:,3), 'Color', c(round(tt/t*256),:));
    end
end
% axis([-100 900 0 1000]);
set(gca, 'Color', 'none');

%% eliminate repeated links
epsilon = 5;
remove_cnt = 0;
for ff = 2:4
    ff
    for tt = 1:t-1
        remove_flag = zeros(size(matchingAll{ff}{tt},1), 1);
        for ii = 1:size(matchingAll{ff}{tt},1)/3
            loc1 = matchingAll{ff}{tt}(ii*3-2,:);
            for pp = 1:ff-1
                for jj = 1:size(matchingAll{pp}{tt},1)/3
                    loc2 = matchingAll{pp}{tt}(jj*3-2,:);
                    if norm(loc1 - loc2) < epsilon
                        if norm(matchingAll{ff}{tt}(ii*3-1,:) - ...
                                matchingAll{pp}{tt}(jj*3-1,:)) < epsilon
                            remove_flag(ii*3-2:ii*3) = 1;
                            remove_cnt = remove_cnt + 1;
                            continue;
                        end
                    end
                end
            end
        end
        matchingAll{ff}{tt}(logical(remove_flag),:) = [];
    end
end