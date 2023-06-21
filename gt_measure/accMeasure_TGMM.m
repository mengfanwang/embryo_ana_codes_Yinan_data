clc;clear;close all;
% Runn TGMM wrapper to get tracking result
gt_path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';

%% get ground truth
% should be consistent with accMeasure.m
link_filename = fullfile(gt_path, 'MastodonTable-Link.csv');
link_table = readtable(link_filename);

link_table = table2array(link_table(3:end, [4 5 7 8 9 10 11 12]));
link_ind = any(link_table(:,3:end),2);
link_table = link_table(link_ind,:);

% plot lineages
G = digraph();
EndNodes = cell(size(link_table,1),2);
for ii = 1:size(link_table,1)
    EndNodes{ii,1} = num2str(link_table(ii,1));
    EndNodes{ii,2} = num2str(link_table(ii,2));
end
G = addedge(G, EndNodes(:,1), EndNodes(:,2));
% plot(G, 'Layout', 'layered');

% get spot location
spot_filename = fullfile(gt_path, 'MastodonTable-Spot.csv');
spot_table = readtable(spot_filename);

spot_table = table2array(spot_table(3:end, [1 4 5 6 7 9 10 11 12 13 14])); % id time x y z
spot_ind = find(any(spot_table(:,6:11),2));
spot_table = spot_table(spot_ind,:);
% convert to matlab corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1 1 1/5.86] + 1;

% check no linked cells
inLink_flag = zeros(size(spot_table,1),1);
for ii = 1:size(spot_table)
    if ~ismember(spot_table(ii,1), link_table(:,1)) && ~ismember(spot_table(ii,1), link_table(:,2))
        inLink_flag(ii) = 1;
    end
end
spot_table(find(inLink_flag),:) = [];

% %% load tgmm result
t = 192;
% % tic;
% svList = cell(t,1);
% addpath('/home/mengfan/ForExecute/embryo_ana_codes_Yinan_data/TGMM_wrapper/readTGMM_XMLoutput');
% file_path = '/work/Mengfan/TGMM2.0/Outputs/Tracking_results/21-01-11/XML_finalResult_lht/GMEMfinalResult_frame';
% for tt = 1:1
%     file_name = [file_path num2str(tt-1,'%.4d') '.svb'];       
%     [svList{tt}, sizeIm] = readListSupervoxelsFromBinaryFile(file_name);
% end

% movieInfo.voxIdx = voxIdx;
%%
movieInfo = struct();

load('/work/Mengfan/TGMM2.0/Outputs/Tracking_results/21-01-11/voxIdx.mat');
movieInfo.voxIdx = voxIdx;

% important feature: 
% 3-5: xyz location (not used actually, use svb)
% 7:   parent
% 8:   time
% 10:  skeleton id
load('/work/Mengfan/TGMM2.0/Outputs/Tracking_results/21-01-11/xml2matlab.mat');
trackingMatrix(:,8) = trackingMatrix(:,8) + 1;

%% trackingMatrix to tracks
% get tracks from trackingMatrix
track_num = max(trackingMatrix(:,10));
tracks = cell(track_num, 1);
for ii = 1:size(trackingMatrix,1)
    parent = trackingMatrix(ii,7);
    track_id = trackingMatrix(ii,10);
    division_flag = false;
    start_flag = false;
    normal_flag = false;
    if parent == -1
        % case 1: start of a track
        if isempty(tracks{track_id})
            tracks{track_id}{1} = ii;
            start_flag = true;
        else
            error('Parent not find!');
        end
    else
        % case 2: a normal child
        for jj = 1:length(tracks{track_id})
            if tracks{track_id}{jj}(end) == parent
                tracks{track_id}{jj}(end+1) = ii;
                normal_flag = true;
            end
        end
    end
    if ~division_flag && ~normal_flag
        % case 3: judge a division 
        for jj = 1:length(tracks{track_id})
            if ismember(parent, tracks{track_id}{jj})
                division_flag = true;
                break;
            end
        end
    end
    if division_flag
        tracks{track_id}{end+1} = [parent ii];
    end
    if ~start_flag && ~normal_flag && ~division_flag
        error('Not definied case!');
    end
end
% reshape tracks
movieInfo.tracks = cell(1e6, 1);
track_num = 0;
for ii = 1:length(tracks)
    for jj = 1:length(tracks{ii})
        track_num = track_num + 1;
        movieInfo.tracks{track_num} = tracks{ii}{jj};
    end
end
movieInfo.tracks = movieInfo.tracks(1:track_num);
% get frame info
movieInfo.frames = trackingMatrix(:,8);
movieInfo.n_perframe = zeros(t,1);
for tt = 1:t
    movieInfo.n_perframe(tt) = sum(trackingMatrix(:,8) == tt);
end

%% measure detection accuracy
data_size = [1920 1920 180];
cell_num = size(spot_table,1);
detection_id = nan(cell_num,1);
detection_frame = zeros(cell_num,1);
fprintf('Measure detetection accuracy...\n');
for ii = 1:cell_num
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, cell_num);
    end
    ii_loc = round(spot_table(ii, 3:5));
    ii_idx = sub2ind(data_size, ii_loc(2), ii_loc(1), ii_loc(3));
    ii_time = spot_table(ii, 2);
    ind_candi = find(movieInfo.frames == ii_time);
    ind_match = find(cellfun(@(x) ismember(ii_idx, x), movieInfo.voxIdx(ind_candi)));
    if isempty(ind_match)
        warning('Missing cell.');
    else
        detection_id(ii) = ind_match + sum(movieInfo.n_perframe(1:ii_time-1));
        detection_frame(ii) = ii_time;
    end
end
fprintf('Missing cell num: %d ratio: %f %%\n', sum(isnan(detection_id)), sum(isnan(detection_id))/cell_num*100);
fprintf('Over-merge num: %d ratio: %f %%\n', (cell_num - length(unique(detection_id))), (cell_num - length(unique(detection_id)))/cell_num*100);

%% measure tracking accuracy
link_num = size(EndNodes,1);
link_degree = outdegree(G);

division_cell = G.Nodes.Name(link_degree==2);
division_cell = cell2mat(cellfun(@str2num, division_cell, 'UniformOutput',false));
for ii = 1:length(division_cell)
    division_cell(ii) = find(spot_table(:,1) == division_cell(ii));
end
division_count = 0;
division_list = zeros(link_num, 12);
missDetect_count = 0;
jump_count = 0;
wrongLink_count = 0;
break_count = 0;
missDetect_list = zeros(link_num, 12);
jump_list = zeros(link_num, 12);
wrongLink_list = zeros(link_num, 12);
break_list = zeros(link_num, 12);
% table info [2 8 2]:
% [s/t index in movieInfo    s/t txyz location    s/t index in mastodon] 
fprintf('Measure tracking accuracy...\n');
tic;
for ii = 1:link_num    
    if mod(ii,1000) == 0
        fprintf('%d / %d\n', ii, link_num);
    end
    s_ind = str2num(G.Edges.EndNodes{ii,1});
    s_ind = find(spot_table(:,1) == s_ind);
    s_corrd = spot_table(s_ind, 2:5);
    t_ind = str2num(G.Edges.EndNodes{ii,2});
    t_ind = find(spot_table(:,1) == t_ind);
    t_corrd = spot_table(t_ind, 2:5);
    if ismember(s_ind, division_cell)
        division_flag = 1;
    else
        division_flag = 0;
    end
    s_ind_old = str2num(G.Edges.EndNodes{ii,1});
    t_ind_old = str2num(G.Edges.EndNodes{ii,2});
    s_ind = detection_id(s_ind);
    t_ind = detection_id(t_ind);
    
    find_flag = 0;
    error_flag = 1;
    missDetect_flag = 0;
    jump_flag = 0;
    wrongLink_flag = 0;
    break_flag = 0;
    if isnan(s_ind) || isnan(t_ind)
        missDetect_flag = 1;
    else
        for jj = 1:length(movieInfo.tracks)
            if ~ismember(s_ind, movieInfo.tracks{jj})
                continue;
            else
                find_flag = 1;
                track = movieInfo.tracks{jj};
                s_loc = find(track == s_ind);
                if length(track) > s_loc
                    if track(s_loc+1) == t_ind
                        error_flag = 0;             % correct
                        break;
                    elseif movieInfo.frames(track(s_loc+1)) > ...
                            movieInfo.frames(track(s_loc)) + 1
                        jump_flag = 1;              % jump
                    else
                        wrongLink_flag = 1;         % wrong link
                    end
                else
                    break_flag = 1;                 % break
                end
            end
        end
        if ~find_flag    % parent is detected but not in any tracks
        for jj = 1:length(movieInfo.tracks)
            if ~ismember(t_ind, movieInfo.tracks{jj})
                continue;
            else
                find_flag = 1;
                track = movieInfo.tracks{jj};
                t_loc = find(track == t_ind);
                if t_loc > 1
                    if movieInfo.frames(track(t_loc-1)) < ...
                            movieInfo.frames(track(t_loc)) - 1
                        jump_flag = 1;              % jump
                    else
                        wrongLink_flag = 1;         % wrong link
                    end
                else
                    break_flag = 1;                 % break
                end
            end
        end           
        end
    end
    if error_flag
%         warning('Error found.');
        if division_flag
            division_count = division_count + 1;
            division_list(division_count, :) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        elseif missDetect_flag
            missDetect_count = missDetect_count + 1;
            missDetect_list(missDetect_count, :) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        elseif jump_flag
            jump_count = jump_count + 1;
            jump_list(jump_count, :) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        elseif wrongLink_flag
            wrongLink_count = wrongLink_count + 1;
            wrongLink_list(wrongLink_count,:) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        elseif break_flag
            break_count = break_count + 1;
            break_list(break_count, :) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        else
            warning('Not considered case found. Added to miss detection.');
            missDetect_count = missDetect_count + 1;
            missDetect_list(missDetect_count, :) = [s_ind t_ind s_corrd t_corrd s_ind_old t_ind_old];
        end
    end
end
toc
division_list = division_list(1:division_count, :);
missDetect_list = missDetect_list(1:missDetect_count, :);
jump_list = jump_list(1:jump_count, :);
wrongLink_list = wrongLink_list(1:wrongLink_count, :);
break_list = break_list(1:break_count, :);
fprintf('Error Summay:\n');
fprintf('Division num:%d ratio:%f %%\n', division_count, division_count/link_num*100);
fprintf('Miss detection num:%d ratio:%f %%\n', missDetect_count, missDetect_count/link_num*100);
fprintf('Jump num:%d ratio:%f %%\n', jump_count, jump_count/link_num*100);
fprintf('Wrong link num:%d ratio:%f %%\n', wrongLink_count, wrongLink_count/link_num*100);
fprintf('Break num:%d ratio:%f %%\n', break_count, break_count/link_num*100);