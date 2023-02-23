function I = display_bnd_map(im, mask, h)
% if nargin < 3
%     fig_on = true;
% end
ids = unique(mask(mask>0));
num_seg = length(ids);
b_all = cell(num_seg, 1);
for i=1:num_seg
    b_all{i} = bwboundaries(mask==ids(i));
end
b_all = cat(1, b_all{:});

% if fig_on
%     h = figure('visible','on');
% else
%     h = figure('visible','off');
% end
if nargin == 2
    h = figure;
end
imshow(im); hold on;
% ha = axes('Parent',h);
for k=1:length(b_all)
   boundary = b_all{k};
   if size(boundary,1) < 10
       continue;
   end
   plot(boundary(:,2), boundary(:,1), 'Color', [0 1 1]/2,'LineWidth',1);
end
hold off;

%I = getimage(h);
I = frame2im(getframe(gcf));
% close(h);
end