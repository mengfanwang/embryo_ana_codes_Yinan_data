function [link_table, spot_table, G] = readGT(path, inLinkOnly)

link_filename = fullfile(path, 'MastodonTable-Link.csv');
link_table = readtable(link_filename);

link_table = table2array(link_table(3:end, [4 5 7 8 9 10 11 12]));
link_ind = any(link_table(:,3:end),2);
link_table = link_table(link_ind,:);
if all(link_table(:,1) > link_table(:,2))
    link_table(:,[1 2]) = link_table(:,[2 1]);
end

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
spot_filename = fullfile(path, 'MastodonTable-Spot.csv');
spot_table = readtable(spot_filename);

% id time x y z
spot_table = table2array(spot_table(3:end, [1 4 5 6 7 9 10 11 12 13 14])); 
spot_ind = find(any(spot_table(:,6:11),2));
spot_table = spot_table(spot_ind,:);

% remove no linked cells
if inLinkOnly
    inLink_flag = zeros(size(spot_table,1),1);
    for ii = 1:size(spot_table)
        if ~ismember(spot_table(ii,1), link_table(:,1)) && ~ismember(spot_table(ii,1), link_table(:,2))
            inLink_flag(ii) = 1;
        end
    end
    spot_table(find(inLink_flag),:) = [];
end

end