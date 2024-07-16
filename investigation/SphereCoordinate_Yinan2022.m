load('tracks_Yinan2022_200_inf_10.mat');

close all;

rx = 840; ry = 840; rz = 1000;
x0 = 960; y0 = 960; z0 = 0;

loc = tracks(tracks(:,8)==0, 3:5);
loc(:,3) = loc(:,3) * 5.86;

% loc = loc - [x0 y0 z0];
% mat = viewmtx(-45,180);
% loc = [loc ones(size(loc,1),1)] * mat;
% loc = loc(:, 1:3);

[X,Y,Z] = sphere;
X = X * rx;
Y = Y * ry;
Z = Z * rz;
figure(1);
plot3(loc(:,1), loc(:,2), loc(:,3), '.'); hold on;
surf(X+x0,Y+y0,Z+z0,'FaceColor','none');

loc_norm = (loc - [x0 y0 z0]) ./ [rx ry rz];
phi = acos(loc_norm(:,3) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2 + loc_norm(:,3).^2));
theta = sign(loc_norm(:,2)) .* acos(loc_norm(:,1) ./ sqrt(loc_norm(:,1).^2 + loc_norm(:,2).^2));
lineage_num = size(loc,1);      
c = colormap('turbo');
lineage_color = c(round(min(phi,2)/2*256),:);

figure(2);
scatter(cos(theta).*phi, sin(theta).*phi, 5, lineage_color, 'filled');

%% cartesein coordinate

    writerObj = VideoWriter('myVideo.avi');
    writerObj.FrameRate = 30;
    open(writerObj);
for tt = 0:190
    tt
ind = find(tracks(:,8)==tt);

lineage = zeros(length(ind),1);
for ii = 1:length(ind)
    if tracks(ind(ii),10) <= lineage_num
        lineage(ii) = tracks(ind(ii),10);
    end
end

loc = tracks(ind, 3:5);
loc = loc(lineage>0,:);
lineage = lineage(lineage>0);

figure(3);
% scatter(-loc(:,2)+2000, loc(:,3)*5.86, 18, lineage_color(lineage,:), 'filled');
h = scatter3(loc(:,1), loc(:,2), loc(:,3)*5.86, 18, lineage_color(lineage,:), 'filled');
pts = [h.XData; h.YData; h.ZData; ones(size(h.XData))];
mat = viewmtx(45,40);
pt2 = mat * pts;

figure(4);
scatter(pt2(1,:)./pt2(4,:), pt2(2,:)./pt2(4,:), 18, lineage_color(lineage,:), 'filled');
axis([300 2200 -600 1100]);
% axis([0 3000 -1000 5000]);
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