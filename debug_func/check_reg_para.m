home_folder = fullfile('/home', getenv('USER'));
addpath(fullfile(home_folder, 'Dropbox/cc_ImHandle/'));


i1 = tifread("/work/public/sameViewFusion/sameViewFusion_050-149_11/00081.tif");
i2 = tifread("/work/public/sameViewFusion/sameViewFusion_050-149_11/00082.tif");



[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 100;
output = cell(1000);
cnt = 1;
for i = 1:200:1400
    for j = 1:200:1400
        for k = 1:20:140
            ii1 = i1(i:i+199, j:j+199, k:k+19);
            ii2 = i2(i:i+199, j:j+199, k:k+19);
            tform = imregtform(ii2,ii1, ...
    "rigid",optimizer,metric);
            disp(tform.T);
            output{cnt} = tform;
            cnt = cnt + 1;
        end
    end
end
