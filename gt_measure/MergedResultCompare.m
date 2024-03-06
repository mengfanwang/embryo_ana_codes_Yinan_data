clc;clear;close all;

%% load data
gt_path = '/work/Nova/embryo_res_folder/mengfan_data_res/Yinan_curated_000_191_0527/';
[link_table, spot_table, G] = readGT(gt_path, 0);
spot_table(logical(spot_table(:,8) & ~spot_table(:,6) & ~spot_table(:,7)& ...
    ~spot_table(:,9)& ~spot_table(:,10)& ~spot_table(:,11)),:) = [];

merged_path = '/work/Nova/embryo_res_folder/mengfan_data_res/gt2result/';
[link_table2, spot_table2, G2] = readGT(merged_path, 0);

spot_table(:, 8) = [];
spot_table2(:, 8) = [];
%% check grond truth result
for ii = 6:10
    a = find(spot_table(:,7));
    a = spot_table(a,1:5);
    b = find(spot_table2(:,7));
    b = spot_table2(b,1:5);

    if size(a,1) ~= size(b,1)
        error('Wrong number');
    else
        dist_pair = pdist2(a(:,3:5), b(:,3:5));
        if max(min(dist_pair)) > 1
            error('Cannot find match cell');
        end
    end
end