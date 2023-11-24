function [tp, fp, fn, tp_list] = evaluation(gt_node, gt_edge, node, edge, resolution)
% evalutae sparse tracking performance with given ground truth using biparite
% matching, consider as the same cell if distance < radius
% INPUT:
% gt_node: ground truth cell, n1*5, [id x y z t]
% gt_edge: ground truth link, m1*2, [id1 id2]
% node: detected cell, n2*5, [id x y z t radius] 
% edge: detected link, n2*5, [id1 id2]
% resolution: 1*3 data resolution

% build graph for bipirate max matching
tt_list = unique(gt_node(:,5));
gt_node_num = size(gt_node,1);
gt_edge_num = size(gt_edge,1);
transition_arcs = zeros(gt_node_num*10,3);
arc_num = 0;
for ii = 1:length(tt_list)
    tt = tt_list(ii);
    gt_node_cur = find(gt_node(:,5) == tt);
    node_cur = find(node(:,5) == tt);
    dist_pair = pdist2(gt_node(gt_node_cur, 2:4), double(node(node_cur, 2:4)).*resolution);
    [loc_x, loc_y] = ind2sub(size(dist_pair), find(dist_pair<node(node_cur,6)'));
    for jj = 1:length(loc_x)
        arc_num = arc_num + 1;
        transition_arcs(arc_num, :) = [gt_node_cur(loc_x(jj)) node_cur(loc_y(jj)) dist_pair(loc_x(jj), loc_y(jj))];
    end
end
transition_arcs = transition_arcs(1:arc_num, :);
transition_arcs(:,2) = transition_arcs(:,2) + gt_node_num;
transition_arcs(:,3) = transition_arcs(:,3) - max(transition_arcs(:,3)) - 1;
total_node_num = size(gt_node,1)+size(node,1);
detection_arcs = [(1:total_node_num)' zeros(total_node_num,1) ...
                   zeros(total_node_num,1) 1e-6*ones(total_node_num,1)];
[trajectories, ~] = mcc4mot(detection_arcs,transition_arcs);
gt2detect = nan(gt_node_num, 1);
for ii = 1:length(trajectories)
    gt2detect(trajectories{ii}(1)) = trajectories{ii}(2) - gt_node_num;
end

% check edge accuracy
if max(gt_node(:,1)) > 1e8
    error('Too large id in ground truth! Please keep id sequential.');
end
gtID2loc = zeros(max(gt_node(:,1)),1);
for ii = 1:length(gt_node)
    gtID2loc(gt_node(ii,1)) = ii;
end
tp = 0; fp = 0; tp_list = zeros(gt_edge_num,1);
% criteria:
% tp: two nodes are matched
% fp: only one nodes are matched
% fn: no nodes are matched, fn = gt_node_num - fp
% other edges don't influence results
for ii = 1:gt_edge_num 
    if mod(ii, 1000) == 0
        fprintf('Validate %d/%d cells\n', ii, gt_edge_num );
    end
    edge_id = gt2detect(gtID2loc(gt_edge(ii,:)));
    if sum(isnan(edge_id)) == 1
        edge_id = node(edge_id(~isnan(edge_id)),1);
        if ismember(edge_id, edge)
            fp = fp + 1;
        end
    elseif sum(isnan(edge_id)) == 0
        edge_id = node(edge_id,1)';
        if ismember(edge_id, edge, "rows")
            tp = tp + 1;
            tp_list(ii) = 1;
        elseif any(ismember(edge_id, edge))
            fp = fp + 1;
        end
    end
end
fn = gt_node_num - tp;
