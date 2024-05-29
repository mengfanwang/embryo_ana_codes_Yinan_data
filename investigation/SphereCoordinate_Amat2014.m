load('tracks_stitch_100_1_10_224.mat');

close all;

loc = tracks(tracks(:,8)==0, 3:5);
loc(:,3) = -loc(:,3);

rx = 780; ry = 780; rz = 90;
x0 = 900; y0 = 920; z0 = -110;
[X,Y,Z] = sphere;
X = X * rx;
Y = Y * ry;
Z = Z * rz;
% figure(1);
% plot3(loc(:,1), loc(:,2), loc(:,3), '.'); hold on;
% surf(X+x0,Y+y0,Z+z0,'FaceColor','none');

loc_norm = (loc - [x0 y0 z0]) ./ [rx ry rz];
phi = acos(loc_norm(:,3) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2 + loc_norm(:,3).^2));
theta = sign(loc_norm(:,2)) .* acos(loc_norm(:,1) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2));
lineage_num = size(loc,1);
c = colormap('turbo');
lineage_color = c(round(min(phi,2)/2*256),:);

figure(2);
scatter(cos(theta).*phi, sin(theta).*phi, 5, lineage_color, 'filled');

%% map to shperical coordinate

%     writerObj = VideoWriter('myVideo.avi');
%     writerObj.FrameRate = 30;
%     open(writerObj);

for tt = 0:0
    tt
ind = find(tracks(:,8)==tt);

lineage = zeros(length(ind),1);
for ii = 1:length(ind)
    if tracks(ind(ii),10) <= lineage_num
        lineage(ii) = tracks(ind(ii),10);
    end
end

loc = tracks(ind, 3:5);
loc(:,3) = -loc(:,3);
loc = loc(lineage>0,:);
loc_norm = (loc - [x0 y0 z0]) ./ [rx ry rz];
phi = acos(loc_norm(:,3) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2 + loc_norm(:,3).^2));
theta = sign(loc_norm(:,2)) .* acos(loc_norm(:,1) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2));
lineage = lineage(lineage>0);

figure(3);
ps = polarscatter(theta, phi/(2*pi)*180, 18, lineage_color(lineage,:), 'filled');
rlim([0 100]);
set(gcf,'position',[100, 100, 960, 960]);

%     frame = getframe(gcf);
%     frame = frame.cdata;
%     if tt == 0
%         sz = size(frame);
%         sz = sz(1:2);
%     else
%         frame = imresize(frame, sz);
%     end
% 
% 
%     writeVideo(writerObj, frame);
%     close all

end
%     close(writerObj);
%% cartesein coordinate

    writerObj = VideoWriter('myVideo.avi');
    writerObj.FrameRate = 30;
    open(writerObj);

for tt = 0:399
    tt
ind = find(tracks(:,8)==tt);

lineage = zeros(length(ind),1);
for ii = 1:length(ind)
    if tracks(ind(ii),10) <= lineage_num
        lineage(ii) = tracks(ind(ii),10);
    end
end

loc = tracks(ind, 3:5);
loc(:,3) = -loc(:,3);
loc = loc(lineage>0,:);
% loc_norm = (loc - [x0 y0 z0]) ./ [rx ry rz];
% phi = acos(loc_norm(:,3) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2 + loc_norm(:,3).^2));
% theta = sign(loc_norm(:,2)) .* acos(loc_norm(:,1) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2));
lineage = lineage(lineage>0);

figure(3);
scatter(loc(:,1), loc(:,3)*5.86+1400, 18, lineage_color(lineage,:), 'filled');
axis([0 1800 0 1400]);
set(gcf,'position',[100, 100, 960, 960]);

    frame = getframe(gcf);
    frame = frame.cdata;
    if tt == 0
        sz = size(frame);
        sz = sz(1:2);
    else
        frame = imresize(frame, sz);
    end


    writeVideo(writerObj, frame);
    close all

end
    close(writerObj);