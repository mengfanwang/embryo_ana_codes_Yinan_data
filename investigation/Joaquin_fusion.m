clc;clear;close all;
addpath ../gt_measure/  ../src_code_matlab/ ../TGMM_wrapper/
addpath /home/mengfan/ForExecute/cc_ImHandle/
dbstop if error

save_path = '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Tracking/Fusion_000_019/';

%% load registraiton info from xml
xml_name = "/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/20220518 isl2b H2Bmcherry overnight.xml";
t = 20;
num_view = 16;
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
% [h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));
% num_time = length(h5_struct);
% num_total = num_time * num_view;

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


%% tif2bdv with registration

% view_name = {'8', '9', '10', '11'};
% timepts_to_process = generate_tps_str(000:019);
% tif_folder_path = {'/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view9', ...
%                    '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view10', ...
%                    '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view11', ...
%                    '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view12'};
% save_data_name = [save_path 'embryo_data_h5'];
% st_loc = []; sz_crop = [];
% tif2bdv_globalReg_multiview(view_name, tif_folder_path, save_data_name, timepts_to_process, st_loc, sz_crop, trans);

%%  mat2tgmm with registration
movieInfo_paths = {'/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Tracking/view9_0218_0_19', ...
                   '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Tracking/view10_0218_0_19', ...
                   '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Tracking/view11_0221_0_19', ...
                   '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Tracking/view12_0214_0_19'};
xml_folder = fullfile(save_path, 'tgmm_mat');
if ~exist(xml_folder)
    mkdir(xml_folder);
end
node_num= 0;
lineage_num = 0;
sc_f = [2, 2, 1];
docNode = cell(t,1);
root = cell(t,1);
z_resolution = trans{1}{1}(3,3);
for vv = 1:4
    load(fullfile(movieInfo_paths{vv}, 'movieInfo.mat'));

    
    start_ind = [0; cumsum(movieInfo.n_perframe)];
    n_total = length(movieInfo.xCoord);
    lineage_list = nan(n_total,1);
    x_ind = (movieInfo.xCoord - 1) * sc_f(1);
    y_ind = (movieInfo.yCoord - 1) * sc_f(2);
    z_ind = (movieInfo.zCoord - 1) * sc_f(3);  
    
    for ii = 1:length(movieInfo.tracks)
        lineage_list(movieInfo.tracks{ii}) = ii - 1;
    end
    for tt = 1:t
        if isempty(docNode{tt})
            docNode{tt} = com.mathworks.xml.XMLUtils.createDocument('document');
            root{tt} = docNode{tt}.getDocumentElement;
        end
        for ii = start_ind(tt)+1:start_ind(tt+1) 
            lineage = lineage_list(ii);
            if isnan(lineage)
                continue;
            end
            id = ii;
            parent = -1;
            if ~isempty(movieInfo.parents{ii})
                parent = movieInfo.parents{ii};
            end
    
            gmm = docNode{tt}.createElement('GaussianMixtureModel');
            gmm.setAttribute('id',num2str(id + node_num));
            gmm.setAttribute('lineage',num2str(lineage + lineage_num));
            gmm.setAttribute('parent',num2str(parent + node_num));
            gmm.setAttribute('splitScore','3');
            gmm.setAttribute('scale','1 1 1');
            gmm.setAttribute('nu','100');
            gmm.setAttribute('beta','100');
            gmm.setAttribute('alpha','100');
            ind_tmp = trans{tt}{vv} * [y_ind(ii) x_ind(ii) z_ind(ii) 1]';
            gmm.setAttribute('m',[num2str(ind_tmp(1)) ' ' num2str(ind_tmp(2)) ' ' num2str(ind_tmp(3)/z_resolution)]);
            gmm.setAttribute('W','0.01 0 0 0 0.01 0 0 0 0.01');
            gmm.setAttribute('svIdx',num2str(id + node_num));
            gmm.appendChild(docNode{tt}.createComment('.'));
            root{tt}.appendChild(gmm);
        end
    end
    node_num= node_num + n_total;
    lineage_num = lineage_num + length(movieInfo.tracks);
end
for tt = 1:t
    ind = char(sprintf("%04d", tt-1));
    xmlwrite(fullfile(xml_folder, ['GMEMfinalResult_frame' ind '.xml']), docNode{tt});
end

