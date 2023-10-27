function detect_splitPair = component_divDetection(movieInfo, paras)
    t = length(movieInfo.n_perframe);
    data_size = paras.data_size;
    cell_num = length(movieInfo.xCoord);
    
    %% find candidates: two components and satisfy size constraint
    div_num = 0;
    detect_splitPair = zeros(cell_num,2);
    data_template = zeros(data_size);
    fprintf('Component detection step 1: start to check components of detections...\n');
    for cell_id = 1:cell_num
        if mod(cell_id,10000) == 0
            fprintf('%d / %d\n', cell_id, cell_num);     
        end
        [vBin, vIdx, ~, ~] = crop3D(data_template, movieInfo.voxIdx{cell_id}, [0 0 0]);
        vBin(ismember(vIdx, movieInfo.voxIdx{cell_id})) = 1;
        cc = bwconncomp(vBin, 6);
        if cc.NumObjects == 2 && ~isempty(movieInfo.parents{cell_id})
            size_ratio = cellfun(@length, cc.PixelIdxList);
            size_ratio = max(size_ratio) / min(size_ratio);
            if size_ratio < paras.size_ratio
                div_num = div_num+1;
                detect_splitPair(div_num,:) = [movieInfo.parents{cell_id} cell_id];
            end
        end
    end
    detect_splitPair = detect_splitPair(1:div_num,:);
    % remove duplicated cases
    detect_splitPair(ismember(detect_splitPair(:,1), detect_splitPair(:,2)),:) = [];
    div_num = size(detect_splitPair,1);

    %% only the head of a track can be child
    head_list = nan(length(movieInfo.tracks),1);
    for ii = 1:length(movieInfo.tracks)
        track = movieInfo.tracks{ii};
        if length(track) >= 1
            head_list(ii) = track(1);
        end
    end
    
    %% filtering: has two head kids and satisfy edge constraint
    has_child_flag = zeros(div_num,1);
    child_pair = nan(div_num,2);
    fprintf('Component detection step 2: candidate filtering...\n');
    for ii = 1:div_num
        cell_id = detect_splitPair(ii,2);
        cell_time = movieInfo.frames(cell_id);
        if cell_time < t
            [neighbors, ~] = divFun.findNeighbor(movieInfo, cell_id, paras);
            neighbors(~ismember(neighbors, head_list)) = []; % keep head only   
            if ~isempty(movieInfo.kids{cell_id}) && ~isempty(neighbors)
                % case1: already have a kids
                child_pair(ii,1) = movieInfo.kids{cell_id};
                child_time = movieInfo.frames(child_pair(ii,1));
                if child_time == cell_time + 1
                    % find an extra kid and test edge ratio
                    child_pair(ii,2) = neighbors(1);
                    has_child_flag(ii) = divFun.edge_ratio_test(movieInfo, cell_id, ...
                        child_pair(ii,:), paras);  
                end
            elseif isempty(movieInfo.kids{cell_id}) && length(neighbors) >= 2
                % case2: no kids
                child_pair(ii,:) = neighbors(1:2);
                has_child_flag(ii) = divFun.edge_ratio_test(movieInfo, cell_id, ...
                        child_pair(ii,:), paras);  
            end
        end
    end
    
    detect_splitPair = [detect_splitPair(:,2) child_pair];
    detect_splitPair(~logical(has_child_flag),:) = [];

end




