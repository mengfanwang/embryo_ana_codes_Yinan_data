clc;clear;close all;

%% get ground truth
addpath ../gt_measure/
path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(path, 1);
% link table: s t
% spot table: id t y x z
link_table = link_table(:,1:2);
spot_table = spot_table(:,1:5);

% convert to tracking result corrdinate
spot_table(:,2) = spot_table(:,2) + 1;
spot_table(:,3:5) = spot_table(:,3:5).*[1/2 1/2 1/5.86] + 1;

%% convert to tracks
d_in = indegree(G);
d_out = outdegree(G);

component_id_list = conncomp(G, 'Type', 'weak')';
tracks = cell(max(component_id_list),1);
link_flag = zeros(size(link_table,1),1);

for component_id = 1:max(component_id_list)
    component_id
    tracks{component_id} = {};
    track_id = 0;
    s_ind = find(component_id_list == component_id & d_in == 0);
    s_node = G.Nodes.Name{s_ind};
    t_ind_list = find(component_id_list == component_id & d_out == 0);
    t_node_list = G.Nodes.Name(t_ind_list);
    for ii = 1:length(t_node_list)
        track_id = track_id + 1;
        tracks{component_id}{track_id} = [];
        path = allpaths(G, s_node, t_node_list{ii});
        if length(path) > 2
            error('Multi paths!');
        end
        path = path{1};
        for jj = 1:length(path) - 1
            parent = str2num(path{jj});
            child = str2num(path{jj+1});
            link_ind = find(ismember(link_table, [parent child], 'rows'));
            if link_flag(link_ind) == 0
                if isempty(tracks{component_id}{track_id})
                    tracks{component_id}{track_id} = [parent child];
                else
                    tracks{component_id}{track_id}(end+1) = child;
                end
                link_flag(link_ind) = 1;
            end
        end
    end
end

%% find division links
ind_in_graph = zeros(max(spot_table(:,1)),1);  % node index in G.Nodes
for ii = 1:length(G.Nodes.Name)
    ind = str2num(G.Nodes.Name{ii});
    ind_in_graph(ind) = ii;
end

division_list = G.Nodes.Name(find(d_out == 2));
division_list(109) = [];   % remove one outlier
div_tracks = cell(length(division_list),3);   % parent child1 child2
for ii = 1:length(division_list)
    parent_id = str2num(division_list{ii});
    component_id = component_id_list(ind_in_graph(parent_id));
    for jj = 1:length(tracks{component_id})
        parent_loc = find(parent_id == tracks{component_id}{jj});
        if ~isempty(parent_loc)
            if parent_loc > 1    % parent track exist
                if isempty(div_tracks{ii,1})
                    div_tracks{ii,1} = tracks{component_id}{jj}...
                        (max(1,parent_loc-9):parent_loc);
                end
            end
            if parent_loc < length(tracks{component_id}{jj}) % child track exist
                child_track = tracks{component_id}{jj}(parent_loc...
                    :min(length(tracks{component_id}{jj}),parent_loc+9));
                if isempty(div_tracks{ii,2})
                    div_tracks{ii,2} = child_track;
                elseif isempty(div_tracks{ii,3})
                    div_tracks{ii,3} = child_track;
                else
                    error('Undefined case!');
                end
            end
        end
    end
end

%% load data
timepts_to_process = generate_tps_str(0:191);
data_folder = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_10';
sc_f = 2;
st_loc = [];
sz_crop = [];

tif_files = dir(fullfile(data_folder, '*.tif'));
if ~isempty(timepts_to_process)
    tif_files(numel(timepts_to_process)+1:end) = [];
    for f = 1:numel(timepts_to_process)
        tif_files(f).name = timepts_to_process(f) + '.tif';
    end     
end

tic;
fprintf('Reading data...');
scale_term = 500;
embryo_vid = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    i
    embryo_vid_temp{1} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid_temp{1} = 255*embryo_vid_temp{1}./scale_term;
    [~, embryo_vid_temp, ~, ~, ~, ~] = data_scaling(sc_f, st_loc, ...
    sz_crop, {}, embryo_vid_temp, {}, {}, {}, {});
    embryo_vid{i} = embryo_vid_temp{1};
end
toc

%% load movieInfo
load('/work/Nova/embryo_res_folder/mengfan_data_res/view10_0504_000_191/movieInfo.mat');

%% find corresponding index
data_size = [960 960 181];
cell_num = size(spot_table,1);
detection_id = cell(size(div_tracks)); 
fprintf('Find cells in movieInfo...\n');
for ii = 1:length(division_list)
    for jj = 1:3
        detection_id{ii,jj} = nan(size(div_tracks{ii,jj}));  
        for kk = 1:length(div_tracks{ii,jj})
            cell_id = find(spot_table(:,1) == div_tracks{ii,jj}(kk));
            cell_loc = round(spot_table(cell_id, 3:5));
            cell_idx = sub2ind(data_size, cell_loc(2), cell_loc(1), cell_loc(3));
            cell_time = spot_table(cell_id, 2);
            ind_candi = find(movieInfo.frames == cell_time);
            ind_match = find(cellfun(@(x) ismember(cell_idx, x), movieInfo.voxIdx(ind_candi)));
            if isempty(ind_match)
                warning('Missing cell.');
            else
                detection_id{ii,jj}(kk) = ind_match + sum(movieInfo.n_perframe(1:cell_time-1));
            end
        end
    end
end

%% get the null size change 
t = 192;
null_size = cell(t,1);
null_size_mean = zeros(t,1);
null_size_std = zeros(t,1);
for tt = 1:t
    null_size{tt} = nan(movieInfo.n_perframe(tt),1);
    for ii = 1:movieInfo.n_perframe(tt)
        cell_id = ii + sum(movieInfo.n_perframe(1:tt-1));
        parent_id = movieInfo.parents{cell_id};
        if isempty(parent_id)
            continue;
        end
        grandpa_id = movieInfo.parents{parent_id};
        if isempty(grandpa_id)
            null_size{tt}(ii) = length(movieInfo.voxIdx{cell_id}) /...
                                length(movieInfo.voxIdx{parent_id});
        else
            null_size{tt}(ii) = length(movieInfo.voxIdx{cell_id})*2 /...
        (length(movieInfo.voxIdx{parent_id})+length(movieInfo.voxIdx{parent_id}));   
        end
%         null_size{tt}(ii) = length(movieInfo.voxIdx{cell_id});
    end
    null_size_mean(tt) = mean(null_size{tt},'omitnan');
    null_size_std(tt) = std(null_size{tt}, 'omitnan');
end

%% investigation - size
normalize_flag = 1;
div_size = cell(size(div_tracks)); 
fprintf('Find cells size...\n');
for ii = 1:length(division_list)
    for jj = 1:3
        div_size{ii,jj} = nan(size(div_tracks{ii,jj}));  
        for kk = 1:length(div_tracks{ii,jj})
            if ~isnan(detection_id{ii,jj}(kk))
                div_size{ii,jj}(kk) = length(movieInfo.voxIdx{detection_id{ii,jj}(kk)});
            end
        end
    end
    if normalize_flag % normalize to parent size
        div_size{ii,2} = div_size{ii,2}/mean(div_size{ii,1}(end-1:end), "omitnan");
        div_size{ii,3} = div_size{ii,3}/mean(div_size{ii,1}(end-1:end), "omitnan");
        div_size{ii,1} = div_size{ii,1}/mean(div_size{ii,1}(end-1:end), "omitnan");
    end
end
% slightly adjustment
adjust_flag = zeros(size(division_list));
for iter = 1:2
for ii = 1:length(division_list)
    if div_size{ii,2}(2) + div_size{ii,3}(2) > 1.5*div_size{ii,2}(1)
        adjust_flag(ii) = 1;
        div_size{ii,1}(end+1) = (div_size{ii,2}(2) + div_size{ii,3}(2))/2;
        div_size{ii,2}(1) = [];
        div_size{ii,3}(1) = [];
    end
end
end


% plot size
figure(1); hold on;
total_size = nan(length(division_list)*3,20);
for ii = 1:length(division_list)
    for jj = 1:3
        if jj == 1 % parent track
            plot(1-length(div_intensity{ii,jj}):0, div_intensity{ii,jj});
            total_size((ii-1)*3+jj, 11-length(div_intensity{ii,jj}):10) = div_intensity{ii,jj};
        else
            plot(0:length(div_intensity{ii,jj})-1, div_intensity{ii,jj});
            total_size((ii-1)*3+jj, 10:length(div_intensity{ii,jj})+9) = div_intensity{ii,jj};
        end
    end
end
mean_size = mean(total_size(:,4:14),'omitnan');
std_size = std(total_size(:,4:14),'omitnan');

%% investigation - intensity
% get the average intensity of each time point
t = 192;
null_intensity = cell(t,1);
null_intensity_mean = zeros(t,1);
null_intensity_std = zeros(t,1);
for tt = 1:t
    null_intensity{tt} = zeros(movieInfo.n_perframe(tt),1);
    for ii = 1:movieInfo.n_perframe(tt)
        cell_id = ii + sum(movieInfo.n_perframe(1:tt-1));
        null_intensity{tt}(ii) = median(embryo_vid{tt}(movieInfo.voxIdx{cell_id}));
    end
    null_intensity_mean(tt) = mean(null_intensity{tt});
    null_intensity_std(tt) = std(null_intensity{tt});
end

%% upper tail truncated gaussian fitting
for tt = 1:t
    data = null_intensity{tt};
    upper_threshold = 2*median(data)-min(data);  % symmetric at median  
    phi_beta = 1 - sum(data > upper_threshold)/length(data);  
    data = data(data <= upper_threshold);
    beta = norminv(phi_beta);
    var_ratio = 1 - beta*normpdf(beta)/phi_beta - (normpdf(beta)/phi_beta)^2;
    null_intensity_std(tt) = sqrt(var(data)/var_ratio);
    null_intensity_mean(tt) = mean(data) + null_intensity_std(tt)*normpdf(beta)/phi_beta;
end

%% get the intensity of divisions
figure(2); hold on
total_intensity = nan(length(division_list)*3,20);
div_intensity = cell(size(div_tracks)); 
combine_score = zeros(length(division_list),1);
fprintf('Find cells intensity...\n');
for ii = 1:length(division_list)
    for jj = 1:3
        for kk = 1:length(div_tracks{ii,jj})
            if ~isnan(detection_id{ii,jj}(kk))
                tt = movieInfo.frames(detection_id{ii,jj}(kk));
                div_intensity{ii,jj}(kk) = median(embryo_vid{tt}(movieInfo.voxIdx{detection_id{ii,jj}(kk)}));
                div_intensity{ii,jj}(kk) = div_intensity{ii,jj}(kk) / null_intensity_mean(tt);
            end
        end
        if jj == 1 % parent track
            plot(1-length(div_intensity{ii,jj}):0, div_intensity{ii,jj}, '.');
            total_intensity((ii-1)*3+jj, 11-length(div_intensity{ii,jj}):10) = div_intensity{ii,jj};
        else
            plot(0:length(div_intensity{ii,jj})-1, div_intensity{ii,jj}, '.');
            total_intensity((ii-1)*3+jj, 10:length(div_intensity{ii,jj})+9) = div_intensity{ii,jj};
        end
    end
    parent_inten = div_intensity{ii,1};
    if length(parent_inten) > 3
        parent_inten = parent_inten(end-2:end);
    end
    child_inten1 = div_intensity{ii,2};
    child_inten2 = div_intensity{ii,3};
    if length(child_inten1) > 3
        child_inten1 = child_inten1(2:3);
    else
        child_inten1 = child_inten1(2:end);
    end
    if length(child_inten2) > 3
        child_inten2 = child_inten2(2:3);
    else
        child_inten2 = child_inten2(2:end);
    end
    parent_inten = parent_inten(~isnan(parent_inten));
    child_inten1 = child_inten1(~isnan(child_inten1));
    child_inten2 = child_inten2(~isnan(child_inten2));
    combine_score(ii) = (sum(parent_inten) + sum(child_inten1) + sum(child_inten2))/...
        sqrt(length(parent_inten) + length(child_inten1) + length(child_inten2));
end
mean_intensity = mean(total_intensity(:,4:14),'omitnan');
std_intensity = std(total_intensity(:,4:14),'omitnan');

%% investigation - similarity between cells
figure(3); hold on
div_similarity = cell(size(div_tracks)); 
total_similarity = nan(length(division_list)*3,20);
for ii = 1:length(division_list)
    for jj = 1:3
        for kk = 1:length(div_tracks{ii,jj})-1
            if ~isnan(detection_id{ii,jj}(kk)) && ~isnan(detection_id{ii,jj}(kk+1))
                tt = movieInfo.frames(detection_id{ii,jj}(kk));
                parent_id = detection_id{ii,jj}(kk);
                child_id = detection_id{ii,jj}(kk+1);
                parent_center = round(movieInfo.orgCoord(parent_id,:));
                child_center = round(movieInfo.orgCoord(child_id,:));
                
                parent_left = parent_center - min(movieInfo.vox{parent_id});
                parent_right = max(movieInfo.vox{parent_id}) - parent_center;
                child_left = child_center - min(movieInfo.vox{child_id});
                child_right = max(movieInfo.vox{child_id}) - child_center;
                crop_left = round(mean([parent_left; child_left]));
                crop_right = round(mean([parent_right; child_right]));
                
                parent_im = embryo_vid{movieInfo.frames(parent_id)}(parent_center(2)-crop_left(2):parent_center(2)+crop_right(2),...
                    parent_center(1)-crop_left(1):parent_center(1)+crop_right(1), parent_center(3)-crop_left(3):parent_center(3)+crop_right(3));
                child_im = embryo_vid{movieInfo.frames(child_id)}(child_center(2)-crop_left(2):child_center(2)+crop_right(2),...
                    child_center(1)-crop_left(1):child_center(1)+crop_right(1), child_center(3)-crop_left(3):child_center(3)+crop_right(3));
                coef = sum(parent_im(:).* child_im(:)) / norm(parent_im(:)) / norm(child_im(:));
                div_similarity{ii,jj}(kk) = coef;
            end
        end
        if jj == 1 % parent track
            plot(1-length(div_similarity{ii,jj}):0, div_similarity{ii,jj}, '.');
            total_similarity((ii-1)*3+jj, 11-length(div_similarity{ii,jj}):10) = div_similarity{ii,jj};
        else
            plot(0:length(div_similarity{ii,jj})-1, div_similarity{ii,jj}, '.');
            total_similarity((ii-1)*3+jj, 10:length(div_similarity{ii,jj})+9) = div_similarity{ii,jj};
        end
    end
end
    