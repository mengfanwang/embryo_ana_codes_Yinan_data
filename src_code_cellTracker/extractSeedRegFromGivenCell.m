function [seed_region, missed_frames] = extractSeedRegFromGivenCell(...
    parent_kid_vec, movieInfo, thresholdMaps, embryo_vid)
% given two footprints of the same cell (or one), get the seed in the 
% jumped frame in between
% NOTE: here we do not consider if the seed is already occupied by an
% existing cell, we will check this in the following steps
p_id = parent_kid_vec(1);
k_id = parent_kid_vec(2);
if isnan(p_id) || isnan(k_id)% starting point
    if ~isnan(k_id)
        missed_frames = movieInfo.frames(k_id)-1;
        idx = movieInfo.voxIdx{k_id};
        thres = thresholdMaps{movieInfo.frames(k_id)}(idx);

        % transfer idx with given non-rigid drift info
        sz = size(embryo_vid{movieInfo.frames(k_id)});
        [idx_y, idx_x, idx_z] = ind2sub_direct(size(embryo_vid{movieInfo.frames(k_id)}), idx);
        frame_shift = getNonRigidDrift([0 0 0], [mean(idx_x) mean(idx_y) mean(idx_z)], ...
            missed_frames, movieInfo.frames(k_id), movieInfo.drift, movieInfo.driftInfo);
        idx_x = round(idx_x + frame_shift(1));
        idx_y = round(idx_y + frame_shift(2));
        idx_z = round(idx_z + frame_shift(3));
        % remove pixels out of the view
        invalid_flag = (idx_x < 1) | (idx_x > sz(2)) | (idx_y < 1) | ...
                       (idx_y > sz(1)) | (idx_z < 1) | (idx_z > sz(3));
        idx_x(invalid_flag) = []; 
        idx_y(invalid_flag) = [];
        idx_z(invalid_flag) = [];
        idx = sub2ind_direct(size(embryo_vid{missed_frames}), idx_y, idx_x, idx_z);
    elseif ~isnan(p_id)
        missed_frames = movieInfo.frames(p_id)+1;
        idx = movieInfo.voxIdx{p_id};
        thres = thresholdMaps{movieInfo.frames(p_id)}(idx);
        
        % transfer idx with given non-rigid drift info
        sz = size(embryo_vid{movieInfo.frames(p_id)});
        [idx_y, idx_x, idx_z] = ind2sub_direct(sz, idx);
        % we don't have info in the missed frame. The drift info is
        % approximate by the parent location
        frame_shift = getNonRigidDrift([0 0 0], [mean(idx_x) mean(idx_y) mean(idx_z)], ...
            movieInfo.frames(p_id), missed_frames, movieInfo.drift, movieInfo.driftInfo);
        idx_x = round(idx_x + frame_shift(1));
        idx_y = round(idx_y + frame_shift(2));
        idx_z = round(idx_z + frame_shift(3));
        % remove pixels out of the view
        invalid_flag = (idx_x < 1) | (idx_x > sz(2)) | (idx_y < 1) | ...
                       (idx_y > sz(1)) | (idx_z < 1) | (idx_z > sz(3));
        idx_x(invalid_flag) = []; 
        idx_y(invalid_flag) = [];
        idx_z(invalid_flag) = [];
        idx = sub2ind_direct(size(embryo_vid{missed_frames}), idx_y, idx_x, idx_z);
    else
        error('at least one element in parent_kid_vec should be valid!');
    end
%     if find(thres==0,1)
%         keyboard;
%     end
    
    if isempty(idx)
        seed_region = [];
        return;
    end
    % new criteria to justify seed region based on z-test with acceleration
    [embryo_vid_cropped, ~, ~, loc_org_xyz]  = crop3D(embryo_vid{missed_frames}, idx, [5 5 3]);
    thresholdMaps_cropped = crop3D(thresholdMaps{missed_frames}, idx, [5 5 3]);
    fg = false(size(embryo_vid_cropped));
    idx_original = idx;
    idx = sub2ind_direct(size(embryo_vid_cropped), idx_y-loc_org_xyz(2) + 1, ...
        idx_x-loc_org_xyz(1) + 1, idx_z-loc_org_xyz(3) + 1);
    fg(idx) = true;
    bg = imdilate(fg, ones(5,5,3));
    bg2 = imdilate(fg, ones(3,3,1));
    bg(bg2) = false;
    bg(thresholdMaps_cropped>0) = false;
    bg_vals = embryo_vid_cropped(bg);
    fg(thresholdMaps_cropped>0) = false;
    vals = embryo_vid_cropped(fg);
    z_score = (mean(vals) - mean(bg_vals))/std(bg_vals)*sqrt(length(vals));
    z_threshold = -norminv(0.01/length(movieInfo.vox)); % FDR control
    if z_score < z_threshold
        seed_region = [];
        return;
    else
        seed_region = idx_original;
    end
    
%     % new criteria to justify seed region based on z-test
%     fg = false(size(embryo_vid{missed_frames}));
%     fg(idx) = true;
%     bg = imdilate(fg, ones(5,5,3));
%     bg2 = imdilate(fg, ones(3,3,1));
%     bg(bg2) = false;
%     bg(thresholdMaps{missed_frames}>0) = false;
%     bg_vals = embryo_vid{missed_frames}(bg);
%     fg(thresholdMaps{missed_frames}>0) = false;
%     vals = embryo_vid{missed_frames}(fg);
%     z_score = (mean(vals) - mean(bg_vals))/std(bg_vals)*sqrt(length(vals));
%     z_threshold = -norminv(0.01/length(movieInfo.vox)); % FDR control
%     if z_score < z_threshold
%         seed_region = [];
%         return;
%     else
%         seed_region = idx;
%     end

%     if isempty(thres(thres>0)) %|| thres == 0
%         seed_region = [];
%         return;
%     else
%         thres = mode(thres(thres>0));
%         vals = embryo_vid{missed_frames}(idx);
%         seed_region = idx(vals>thres);
%         if isempty(seed_region)
%             return;
%         end
%     end
else
    missed_frames = movieInfo.frames(p_id)+1 : ...
        movieInfo.frames(k_id)-1;
    % find seed region: can be more complicated
    seed_region = intersect(movieInfo.voxIdx{p_id}, movieInfo.voxIdx{k_id});
    if isempty(seed_region)
        return;
    end
end

%% pick the largest region indicated by the seed region
% way 1
%     tic;
%     bw = false(size(thresholdMaps{missed_frames}));
%     bw(seed_region) = true;
%     out_bw = pickLargestReg(bw, 6);
%     seed_region1 = find(out_bw);
%     toc;
% way 2
[bw, st_pt, idx] = mapIdx2Bw(seed_region, size(thresholdMaps{1}));
[bw, rm_flag, idx] = pickLargestReg(bw, 6, idx);
if rm_flag % otherwise, no need to change seed_region
    [yy, xx, zz] = ind2sub_direct(size(bw), idx);
    seed_region = sub2ind_direct(size(thresholdMaps{1}), yy+(st_pt(1)-1), ...
        xx+(st_pt(2)-1), zz+(st_pt(3)-1));
end
end