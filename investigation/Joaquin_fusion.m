clc;clear;close all;
addpath ../gt_measure/

xml_name = "/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/20220518 isl2b H2Bmcherry overnight.xml";


xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
% [h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));
% num_time = length(h5_struct);
% num_total = num_time * num_view;

num_view = 16;
t = 20;

trans = cell(t, 1);
for tt = 1:t
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
    trans{tt} = cell(4,1);
    for vv = 1:4
        ind = (tt-1)*num_view + (vv+8);
        trans{tt}{vv} = eye(4,4);
        for jj = 1:length(register_info{ind}.ViewTransform)
            trans_temp = cellfun(@str2num, split(register_info{ind}.ViewTransform{jj}.affine.Text));
            trans_temp = [reshape(trans_temp,4,3)'; 0 0 0 1];
            trans{tt}{vv} = trans{tt}{vv}*trans_temp;
        end
%         trans{vv}(1:3,:) = trans{vv}(1:3,:) / downsample_scale;  % downsample
    end
end