function [output1, output2] = matfiles2cell(folder, foldertype, names)
% read the mat files under a given folder. concatnate them into cells.
% if names is given, then only read those mat files with given name.

if nargin == 1
    names = [];
end

files = dir(fullfile(folder, '*.mat'));
if ~isempty(names)
    files(numel(names)+1:end) = [];
    for f = 1:numel(names)
        files(f).name = names(f) + '.mat';
    end
end
output1 = cell(numel(files), 1);
output2 = cell(numel(files), 1);
for f = 1:numel(files)
    if strcmp(foldertype, 'synQuant_refine_res')
        m = load(fullfile(files(f).folder, files(f).name),'refine_res', 'threshold_res');
        output1{f} = m.refine_res;
        output2{f} = m.threshold_res;
    elseif strcmp(foldertype,'synQuant_priCvt_res')
        m = load(fullfile(files(f).folder, files(f).name), 'eig_res_2d', 'eig_res_3d');
        output1{f} = m.eig_res_2d;
        output2{f} = m.eig_res_3d;
    elseif strcmp(foldertype, 'varianceMap')
        m = load(fullfile(files(f).folder, files(f).name),'varMap');
        output1{f} = m.varMap;
    else
        error("undefined file type");
    end
end
end