clc;clear;close all;
% mapping cells in the first time point to the next iteratively
% until the target tp

%% load initial position
load('/work/public/Embryo/Amat2014/Detection/Wei_refine_res/00000.mat');
start_tp = 0;
end_tp = 399;


im_sz = size(refine_res);
cell_num = max(refine_res(:));
cell_loc = zeros(cell_num,3);
for ii = 1:cell_num
    if mod(ii, 100) == 0
        fprintf('%d/%d\n', ii, cell_num);
    end
    cell_idx = find(refine_res == ii);
    [y, x, z] = ind2sub_direct(im_sz, cell_idx);
    cell_loc(ii,:) = [mean(x) mean(y) mean(z)];
end
cell_loc(:,1) = cell_loc(:,1) / 2;
cell_loc(:,2) = cell_loc(:,2) / 2;
boundary_label = cell_loc(:,3) >= 115;

%% plot initial positions
figure(1);title('T = 0');hold on;
scatter3(cell_loc(~boundary_label,1), cell_loc(~boundary_label,2), cell_loc(~boundary_label,3).*8, 36, [0.5 0.5 0.5], '.');
scatter3(cell_loc(boundary_label,1), cell_loc(boundary_label,2), cell_loc(boundary_label,3).*8, 72, [0 0.4470 0.7410], '.');

%% load registration info
timepts_to_process = generate_tps_str(start_tp+1:end_tp);

st_loc = [];
sz_crop = [];
translation_path = '/work/public/Embryo/Amat2014/Registration_reverse';
driftInfo.grid_size = 32;    % vector numbers each dim, currently cube only
driftInfo.batch_size = [30 30 8];

[y_batch, x_batch, z_batch] = meshgrid(0:driftInfo.grid_size+1);
if ~isempty(st_loc)  % adjust if crop
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2 - st_loc(2)/sc_f;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2 - st_loc(1)/sc_f;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2 - st_loc(3);
else
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2;
end
driftInfo.y_batch = y_batch;
driftInfo.x_batch = x_batch;
driftInfo.z_batch = z_batch;  % vector locations in original image with one padding

grid_size = driftInfo.grid_size;
drift.x_grid = cell(length(timepts_to_process),1);
drift.y_grid = cell(length(timepts_to_process),1);
drift.z_grid = cell(length(timepts_to_process),1);
for ii = 1:length(timepts_to_process)
    load(fullfile(translation_path, timepts_to_process(ii)+'.mat'), 'phi_current_vec');
    x_grid = padarray(reshape(phi_current_vec(1:3:end-2), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
    y_grid = padarray(reshape(phi_current_vec(2:3:end-1), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
    z_grid = padarray(reshape(phi_current_vec(3:3:end), [grid_size grid_size grid_size]), [1 1 1], 'replicate'); 
    drift.x_grid{ii} = x_grid;
    drift.y_grid{ii} = y_grid;
    drift.z_grid{ii} = z_grid;
end

%% load mapping result
cell_loc_new = zeros(size(cell_loc));
for ii = 1:size(cell_loc_new,1)
    if mod(ii, 100) == 0
        fprintf('%d/%d\n', ii, cell_num);
    end
    tempVox = cell_loc(ii,:);
    for tt = start_tp+1:200
        y_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.y_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        x_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.x_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        z_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.z_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        tempVox = tempVox + [y_bias x_bias z_bias];
    end
    cell_loc_new(ii,:) = tempVox;
end

%% plot finial positions
figure(2);title('T = 200');hold on;
scatter3(cell_loc_new(~boundary_label,1), cell_loc_new(~boundary_label,2), cell_loc_new(~boundary_label,3).*8, 36, [0.5 0.5 0.5], '.');
scatter3(cell_loc_new(boundary_label,1), cell_loc_new(boundary_label,2), cell_loc_new(boundary_label,3).*8, 72, [1 0 0], '.');

%% load mapping result
cell_loc = cell_loc_new;
cell_loc_new = zeros(size(cell_loc));
for ii = 1:size(cell_loc_new,1)
    if mod(ii, 100) == 0
        fprintf('%d/%d\n', ii, cell_num);
    end
    tempVox = cell_loc(ii,:);
    for tt = 201:end_tp
        y_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.y_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        x_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.x_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        z_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.z_grid{tt}, tempVox(1), tempVox(2), tempVox(3), 'linear');
        tempVox = tempVox + [y_bias x_bias z_bias];
    end
    cell_loc_new(ii,:) = tempVox;
end

%%
figure(3);title('T = 399');hold on;
scatter3(cell_loc_new(~boundary_label,1), cell_loc_new(~boundary_label,2), cell_loc_new(~boundary_label,3).*8, 36, [0.5 0.5 0.5], '.');
scatter3(cell_loc_new(boundary_label,1), cell_loc_new(boundary_label,2), cell_loc_new(boundary_label,3).*8, 72, [1 0 0], '.');
