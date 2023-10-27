function detect_divPair = rule_based_divDetection(movieInfo, embryo_vid, paras)
    t = length(movieInfo.n_perframe);
    data_size = paras.data_size;
    im_resolution = paras.im_resolution;
    length_thre = paras.length_thre;
    child2parent_ratio = paras.child2parent_ratio;
    
    % should related to resolution
    strel_rad = 3;
    strel_ker = ones(strel_rad*2+1, strel_rad*2+1,ceil(strel_rad/3)*2+1);
    [xx,yy,zz] = ind2sub_direct(size(strel_ker), find(strel_ker));
    dist = sqrt( (xx - strel_rad-1).^2 + (yy - strel_rad-1).^2 + ( (zz - ceil(strel_rad/3)-1)*5 ).^2 );
    strel_ker(dist>=strel_rad) = 0;
    strel_ker = strel(strel_ker);

    %% find all detections and candidate divisions in movieInfo
    cell_num = length(movieInfo.xCoord);
    div_candiList = cell(cell_num,1);
    div_candiList_2nd = cell(cell_num,1);
    fprintf('Get head of tracks...\n'); 
    % get heads of tracks
    head_list = nan(length(movieInfo.tracks),1);
    head_list_2nd = nan(length(movieInfo.tracks),1);  
    for ii = 1:length(movieInfo.tracks)
        track = movieInfo.tracks{ii};
        if length(track) >= length_thre
            head_list(ii) = track(1);
            z_inten = zeros(length_thre,1);
            mean_inten = zeros(length_thre,1);
            for jj = 1:length_thre
                voxInten = embryo_vid{movieInfo.frames(track(jj))}(movieInfo.voxIdx{track(jj)});
                mean_inten(jj) = mean(voxInten);
                [vidOut, vIdx, ~, ~] = crop3D(embryo_vid{movieInfo.frames(track(jj))},...
                    movieInfo.voxIdx{track(jj)}, [4 4 2]);
                vBin = zeros(size(vidOut));
                vBin(ismember(vIdx, movieInfo.voxIdx{track(jj)})) = 1;
            
                z_inten(jj) = divFun.bgZtest(vidOut, vBin, strel_ker);
            end
            if z_inten(1) == min(z_inten) && mean_inten(1) == min(mean_inten)
                head_list_2nd(ii) = track(2);
            end
        end
    end
    fprintf('Rule-based detection step 1: find cells and candidate divisions in movieInfo...\n'); 
    tic;
    for parent_id = 1:cell_num
        if mod(parent_id,10000) == 0
            fprintf('%d / %d\n', parent_id, cell_num);
        end
        parent_time = movieInfo.frames(parent_id);
        % if have child already, it should be small
        child_small_flag = true;
        child_id = movieInfo.kids{parent_id};
        if ~isempty(child_id)
            if length(movieInfo.voxIdx{child_id}) >= ...
                length(movieInfo.voxIdx{parent_id}) * child2parent_ratio 
                child_small_flag = false;
            end
        end
        if parent_time ~= t && child_small_flag
            % find the 5-nearest neighbor
            [div_candi, ~] = divFun.findNeighbor(movieInfo, parent_id, paras);
            % only heads are included   
            div_candi1 = div_candi(ismember(div_candi,head_list));
            div_candiList{parent_id} = div_candi1;
    
            % second head included
            div_candi2 = div_candi(ismember(div_candi,head_list_2nd));
            div_candiList_2nd{parent_id} = div_candi2;
        end
    end
    toc                                    
    
    %% convert candiList to division cases
    div_num = 0;
    detect_divPair = zeros(cell_num,4); %[parent child1 child2 new_child(1/2)]
    div_num2nd = 0;
    detect_divPair2nd = zeros(cell_num,4);  
    fprintf('Rule-based detection step 2: filter candidate divisions based on rules...\n'); 
    for ii = 1:cell_num
        if mod(ii,10000) == 0
            fprintf('%d / %d\n', ii, cell_num);
        end
        parent_id = ii;
        detect_flag = false;
        if ~isempty(movieInfo.kids{parent_id}) && ...
                ~isempty(div_candiList{ii})         % case1: one child exists
            child1_id = movieInfo.kids{parent_id};
            children_pair = [repmat(child1_id, [length(div_candiList{ii}) 1])...
                                 div_candiList{ii}'];
            detect_divPair_tmp = divFun.candidateFiltering(movieInfo, parent_id, ...
                                    children_pair, data_size, paras);  
            if ~isempty(detect_divPair_tmp)
                div_num = div_num + 1;
                detect_divPair(div_num,:) = [detect_divPair_tmp 1];
                detect_flag = true;
            end
        elseif length(div_candiList{ii}) >= 2       % case2: no child exists
            children_pair = nchoosek(div_candiList{ii}, 2);
            
            detect_divPair_tmp = divFun.candidateFiltering(movieInfo, parent_id, ...
                                    children_pair, data_size, paras);   
            if ~isempty(detect_divPair_tmp)
                div_num = div_num + 1;
                detect_divPair(div_num,:) = [detect_divPair_tmp 2];
                detect_flag = true;
            end
        end
        if ~detect_flag
            % try to find the second head
            if ~isempty(movieInfo.kids{parent_id})  && ...
                    ~isempty(div_candiList_2nd{parent_id})
                % repeat the same step
                child1_id = movieInfo.kids{parent_id};
                children_pair = [repmat(child1_id, [length(div_candiList_2nd{parent_id}) 1])...
                                    div_candiList_2nd{parent_id}'];
                detect_divPair_tmp = divFun.candidateFiltering(movieInfo, parent_id, ...
                                    children_pair, data_size, paras); 
                if ~isempty(detect_divPair_tmp)
                    div_num2nd = div_num2nd + 1;
                    detect_divPair2nd(div_num2nd,:) = [detect_divPair_tmp 3];
                end
            end
        end
    end
    detect_divPair = detect_divPair(1:div_num,:);
    detect_divPair = unique(detect_divPair, 'rows');
    detect_divPair2nd = detect_divPair2nd(1:div_num2nd,:);
    detect_divPair2nd = unique(detect_divPair2nd, 'rows');
    
    
    %% double check the second head and remove if conflict with the first
    fprintf('Rule-based detection step 3: remove conflict candidates...\n'); 
    new_track_list = zeros(size(detect_divPair,1)*2,1);
    new_track_num = 0;
    for ii = 1:size(detect_divPair,1)
        if detect_divPair(ii,4) == 2
            new_track_num = new_track_num + 1;
            new_track_list(new_track_num) = movieInfo.particle2track(detect_divPair(ii,2),1);
        end
        new_track_num = new_track_num + 1;
        new_track_list(new_track_num) = movieInfo.particle2track(detect_divPair(ii,3),1);
    end
    new_track_list = new_track_list(1:new_track_num);
    remove_flag = zeros(size(detect_divPair2nd,1),1);
    for ii = 1:size(detect_divPair2nd,1)
        if ismember(movieInfo.particle2track(detect_divPair2nd(ii,3),1), new_track_list)
            remove_flag(ii) = 1;
        end
    end
    detect_divPair2nd(logical(remove_flag), :) = [];
    detect_divPair = [detect_divPair; detect_divPair2nd];
    
    
    %% for each pair of the children, mapping their center to the nearest parent
    fprintf('Rule-based detection step 4: futher remove conflicts...\n'); 
    multi_child_list = tabulate(reshape(detect_divPair(:,2:3), [] ,1));
    multi_child_list = find(multi_child_list(:,2) > 1);
    remove_flag = zeros(size(detect_divPair,1),1);
    for ii = 1:size(multi_child_list,1)
        parent_loc = find(any(ismember(detect_divPair(:,2:3), multi_child_list(ii)),2));
        dist2parent = inf(size(parent_loc));
        child_time = movieInfo.frames(multi_child_list(ii));
        for jj = 1:length(parent_loc)
            parent_id = detect_divPair(parent_loc(jj),1);
            child1_id = detect_divPair(parent_loc(jj),2);
            child2_id = detect_divPair(parent_loc(jj),3);
            frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
                movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
            dist2parent(jj) = norm( (mean(movieInfo.vox{parent_id}) + frame_shift -...
                mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}])) .* im_resolution);
        end
        [~, minIdx] = min(dist2parent);
        
        remove_flag(parent_loc) = 1;
        remove_flag(parent_loc(minIdx)) = 0;    
    
    end
    detect_divPair(logical(remove_flag), :) = [];
    % if multiple divisions connect, only keep the first one
    remove_flag = ismember(detect_divPair(:,1),detect_divPair(:,2))...
        | ismember(detect_divPair(:,1),detect_divPair(:,3));
    detect_divPair(logical(remove_flag), :) = [];

end

