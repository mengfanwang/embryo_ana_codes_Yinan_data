
% close all
% load('tracks_stitch_100_1_10.mat');
% tracks(:, 3:5) = tracks(:, 3:5) .* [1 1 8];
% load('tracks_stitch2.mat');
% tracks(:, 3:5) = tracks(:, 3:5) .* [1 1 8];

xlims = [min(tracks(:, 3)), max(tracks(:, 3))];
ylims = [min(tracks(:, 4)), max(tracks(:, 4))];
zlims = [min(-tracks(:, 5)), max(-tracks(:, 5))];   

% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 30;
% open(writerObj);
for tt = 399:399
    tt
    cells = (tracks(:, 8) == tt) & (tracks(:, 9) == 1);
    cellsbk = (tracks(:, 8) == tt) & (tracks(:, 9) == 0);
    cellnew = (tracks(:, 8) == tt) & (tracks(:, 9) == 2);

    % Create figure without displaying itr
%     fig = figure('visible', 'off');

    hold on
    scatter3(tracks(cellsbk, 3), tracks(cellsbk, 4), -tracks(cellsbk, 5),  1, 'o', 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [0.5, 0.5, 0.5]);
    scatter3(tracks(cells, 3), tracks(cells, 4), -tracks(cells, 5), 1, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 3);
    scatter3(tracks(cellnew, 3), tracks(cellnew, 4), -tracks(cellnew, 5), 1, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'b', 'LineWidth', 3);
    axis equal

    view(5, -3); % Change viewing angle
    title(['t = ', num2str(tt)]);

    xlim(xlims);
    ylim(ylims);
    zlim(zlims);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    frame = getframe(gcf);
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
% close(writerObj);