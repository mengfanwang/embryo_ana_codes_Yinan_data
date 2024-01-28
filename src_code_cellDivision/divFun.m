classdef divFun
    % common functions for division detection
    methods(Static)
        function detect_divPair = candidateFiltering(movieInfo, parent_id, ...
                                    children_pair, data_size, paras)  
        % should statisify the following creteria:
        % 0.parent must be the nearest one to the center of the pair of children
        % 1.two child must be the next time points
        % 2.children must split
        % 3.two child size ratio < 3
        % 4.two edge length ratio < 5
            
            im_resolution = paras.im_resolution;
            child_time = movieInfo.frames(parent_id)+1;
            dist2parent = inf(size(children_pair,1),1);
            for jj = 1:size(children_pair,1)
                child1_id = children_pair(jj,1);
                child2_id = children_pair(jj,2);
                frame_shift = getNonRigidDrift([0 0 0], [movieInfo.orgCoord(child1_id,:); ...
                    movieInfo.orgCoord(child2_id,:)], child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
                dist2parent(jj) = norm( (mean(movieInfo.vox{parent_id}) + frame_shift -...
                    mean([movieInfo.vox{child1_id}; movieInfo.vox{child2_id}])) .* im_resolution);
            end
            [~, minIdx] = min(dist2parent);
            child1_id = children_pair(minIdx, 1);
            child2_id = children_pair(minIdx, 2);
        
            adjacent_flag = false;
            if movieInfo.frames(parent_id) + 1 == movieInfo.frames(child1_id) &&...
                movieInfo.frames(parent_id) + 1 == movieInfo.frames(child2_id) 
                adjacent_flag = true;
            end
        
            gap_flag = divFun.division_gap_test(movieInfo, zeros(data_size), child1_id, child2_id);
        
            child_size = [length(movieInfo.voxIdx{child1_id})...
                          length(movieInfo.voxIdx{child2_id})];
            size_flag = (max(child_size) / min(child_size)) < 3;
        
            edge_ratio_flag = divFun.edge_ratio_test(movieInfo, parent_id, ...
                                    [child1_id child2_id], paras);
            
        
            if adjacent_flag && gap_flag && size_flag && edge_ratio_flag
                detect_divPair = [parent_id children_pair(minIdx,:)];
            else
                detect_divPair = [];
            end
        
        end

        function edge_ratio_flag = edge_ratio_test(movieInfo, parent_id, ...
                                    children_pair, paras)  
            im_resolution = paras.im_resolution;
            edge_ratio_constriant = paras.edge_ratio;
            child1_id = children_pair(1);
            child2_id = children_pair(2);
            child_time = movieInfo.frames(child1_id);
            frame_shift1 = getNonRigidDrift([0 0 0], movieInfo.orgCoord(child1_id,:),...
                child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
            euc_dist1 = norm( (movieInfo.orgCoord(parent_id,:) + frame_shift1 -...
                movieInfo.orgCoord(child1_id,:)) .* im_resolution);
            frame_shift2 = getNonRigidDrift([0 0 0], movieInfo.orgCoord(child2_id,:),...
                child_time-1, child_time, movieInfo.drift, movieInfo.driftInfo);
            euc_dist2 = norm( (movieInfo.orgCoord(parent_id,:) + frame_shift2 -...
                movieInfo.orgCoord(child2_id,:)) .* im_resolution);
            edge_ratio = max(euc_dist1, euc_dist2) / min(euc_dist1, euc_dist2);
            edge_ratio_flag = edge_ratio < edge_ratio_constriant;
        end

        function gap_flag = division_gap_test(movieInfo, embryo_vid, child1, child2)
            % get local region
            [vidOut, vIdx, ~, ~] = crop3D(embryo_vid, [movieInfo.voxIdx{child1}; movieInfo.voxIdx{child2}], [0 0 0]);
            vBin = zeros(size(vidOut));
            vBin(ismember(vIdx, movieInfo.voxIdx{child1})) = 1;
            vBin(ismember(vIdx, movieInfo.voxIdx{child2})) = 1;
        
            cc = bwconncomp(vBin, 6);
            if cc.NumObjects == 2
                gap_flag = true;
            else
                gap_flag = false;
            end
        end


        function [neighbors, dist2nei] = findNeighbor(movieInfo, cur_i, paras)
            % find the nearest neighbors in the next frame with given conditions  
            im_resolution = paras.im_resolution;
            max_dist = paras.max_dist;
            max_nei = paras.max_nei;
            
            tt = movieInfo.frames(cur_i);
            curCentroid = movieInfo.orgCoord(cur_i,:);
            nextCentroid = movieInfo.orgCoord(movieInfo.frames==(tt+1),:);
            
            % candidate neighbors roughly selection
            drift = getNonRigidDrift([0 0 0], movieInfo.orgCoord(cur_i,:), tt, tt+1, movieInfo.drift, movieInfo.driftInfo);
            dist_pair = pdist2((curCentroid + drift).*im_resolution, nextCentroid.*im_resolution); % not accurate but enough
            [neighbor_dist, neighbor_candidate] = sort(dist_pair);
            neighbor_dist = neighbor_dist(1,1:min(max_nei, length(dist_pair)));
            neighbor_candidate = neighbor_candidate(1,1:min(max_nei, length(dist_pair)));
            neighbors = neighbor_candidate(neighbor_dist < max_dist);
            neighbors = neighbors + sum(movieInfo.n_perframe(1:tt));
            dist2nei = neighbor_dist(neighbor_dist < max_dist);
        end

        function z_score = bgZtest(data, fMap, strel_ker)
            fg = data(logical(fMap));
            fMap = imdilate(fMap,strel_ker) - fMap;
            bg = data(logical(fMap));
            [mu, sigma] = ordStatApproxKsecWith0s_mat(fg, bg, []);
            L = mean(fg) - mean(bg);
            z_score = (L*sqrt(length(fg)+length(bg)) - mu) / sigma;
        end

        function movieInfo = mergeParent2ChildTrack(movieInfo, parent_id, child_id, child_track)
            child_loc = find(movieInfo.tracks{child_track} == child_id);
            if child_loc == 1
                movieInfo.parents{child_id} = parent_id;
                movieInfo.tracks{child_track} = ...
                    [parent_id; movieInfo.tracks{child_track}];
            elseif child_loc == 2
                movieInfo.parents{child_id} = parent_id;
                movieInfo.tracks{child_track} = ...
                    [parent_id; movieInfo.tracks{child_track}(2:end)];
            else
                error('Wrong case');
            end
        end

    end
end