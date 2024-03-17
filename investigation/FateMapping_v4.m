clc;clear;close all;
% stitch tracking segments
addpath ../src_code_matlab/ ../TGMM_wrapper/

%%
files = {{'000', '049'}, {'040', '089'}};
tps = str2num(files{1}{1}):str2num(files{end}{2});
movieInfoAll.n_perframe = zeros(length(tps),1);
pre_tps = 0;
max_dist = 50;

file_path = '/work/Mengfan/Embryo/Amat2014/Tracking';
translation_path = '/work/public/Embryo/Amat2014/Registration_reverse2';
[driftInfo, drift] = loadRegInfo(translation_path, tps(2:end));

%%
ss = 1;
m1 = load(fullfile(file_path, [files{ss}{1} '_' files{ss}{2}], 'movieInfo.mat'));
m1 = m1.movieInfo;
m2 = load(fullfile(file_path, [files{ss+1}{1} '_' files{ss+1}{2}], 'movieInfo.mat'));
m2 = m2.movieInfo;

%% stitch  two parts: check division and stitch

if pre_tps > tps(1)
    % stitch
end

% try the best to find parents for every cell
% 
cnt = 0;
for ii = 1:length(movieInfo.tracks)
    if isempty(movieInfo.tracks{ii})
        continue
    end
    track_head = movieInfo.tracks{ii}(1);
    if movieInfo.frames(track_head) == 1
        continue
    end
    cur_tps = movieInfo.frames(track_head) - pre_tps;
    cell_loc = movieInfo.orgCoord(track_head);

end

function cell_loc_new = mapCellLoc(cell_loc, start_tp, end_tp, driftInfo, drift)
    cell_loc_new = cell_loc;
    for tt = start_tp:end_tp-1
        y_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.y_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        x_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.x_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        z_bias = interp3(driftInfo.y_batch, driftInfo.x_batch, driftInfo.z_batch,...
            drift.z_grid{tt}, cell_loc_new(:,1), cell_loc_new(:,2), cell_loc_new(:,3), 'linear');
        cell_loc_new = cell_loc_new + [y_bias x_bias z_bias];
    end
end

function [driftInfo, drift] = loadRegInfo(translation_path, tps)
    
    timepts_to_process = generate_tps_str(tps);
    driftInfo.grid_size = 32;    % vector numbers each dim, currently cube only
    driftInfo.batch_size = [30 30 8];
    
    [y_batch, x_batch, z_batch] = meshgrid(0:driftInfo.grid_size+1);
    y_batch = y_batch*driftInfo.batch_size(2) + 0.5 - driftInfo.batch_size(2)/2;
    x_batch = x_batch*driftInfo.batch_size(1) + 0.5 - driftInfo.batch_size(1)/2;
    z_batch = z_batch*driftInfo.batch_size(3) + 0.5 - driftInfo.batch_size(3)/2;
    
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
end