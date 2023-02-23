function out = display_bnd_map_3d(vid, mask)

[h, w, z] = size(vid);
vid = scale_image(vid,0,1);
out1 = cell(z, 1);
h = figure;
for i = 1:z
    disp(i);
    I = display_bnd_map(vid(:,:,i), mask(:,:,i), h);
    out1{i}= I;
end
% out = cat(4, out1{:});
close(h);
minh = inf;
minw = inf;
for i = 1:z
    [hh, w, ~] = size(out1{i});
    minh = min(hh, minh);
    minw = min(w, minw);
end
out = uint8(zeros(minh, minw, 3, z));
for i = 1:z
    [hh, w, ~] = size(out1{i});
    out(:,:,:, i) = out1{i}(hh-minh+1:hh, w-minw+1:w, :);
end

if nargout > 0
    close(h);
else
    zzshow(out);
end
end