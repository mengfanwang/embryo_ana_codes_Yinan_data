function tif2bdv_2channel(tif_folder_path1, tif_folder_path2, save_data_name, timepts_to_process, st_loc, sz_crop)
% tif files to mastodon accepted data format

% clc;clear;close all;
% if isunix
%     addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
%     tif_folder_path = '/work/Mengfan/Embryo/TM0-49/data_reg_0.25';
%     save_data_name = '/work/Mengfan/Embryo/TM0-49/result_0.25';
% else
%     tif_folder_path = 'E:\Embryo\registration_temporal_data\PointReg';
% end

%%% scaling setting %%%
if ~isempty(st_loc)
    end_loc = st_loc + sz_crop - 1;
end

%%% parameter %%%
view_num = 2;     % number of channels
res_mat = [1 1 1; 2 2 2;]';
sub_mat = [16 16 16; 8 8 8]';

% write hdf5 data
h5create([save_data_name '.h5'],'/s00/resolutions',size(res_mat));
h5write([save_data_name '.h5'], '/s00/resolutions', res_mat);
h5create([save_data_name '.h5'],'/s00/subdivisions',size(sub_mat));
h5write([save_data_name '.h5'], '/s00/subdivisions', sub_mat);

h5create([save_data_name '.h5'],'/s01/resolutions',size(res_mat));
h5write([save_data_name '.h5'], '/s01/resolutions', res_mat);
h5create([save_data_name '.h5'],'/s01/subdivisions',size(sub_mat));
h5write([save_data_name '.h5'], '/s01/subdivisions', sub_mat);
l_num = size(res_mat,2);

% read tif file
% tif_files = dir(fullfile(tif_folder_path, '*.tif'));
% time_num = numel(tif_files);
time_num = length(timepts_to_process);
for tt = 1:time_num
    fprintf('Writing time %d\n', tt);
    
    for ss = 0:1
        if ss == 0
            s_ind = '00';
            tif_folder_path = tif_folder_path1;
        elseif ss == 1
            s_ind = '01';
            tif_folder_path = tif_folder_path2;
        end
            
        data = tifread(fullfile(tif_folder_path, timepts_to_process(tt)+'.tif'));
        if ~isempty(st_loc)
            data = data(st_loc(1):end_loc(1), ...
                   st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        end
        [x,y,z] = size(data);
        % data transform
        data_t = zeros(y,x,z);
        for zz = 1:z
            data_t(:,:,zz) = data(:,:,zz)';
        end
        data = data_t;
        data = single(data);
    
        % write h5 file
        t_ind = num2str(100000+tt-1);
        t_ind = t_ind(2:end);
    
        for ll = 0:l_num-1
            l_ind = num2str(ll);
            data_temp = imresize3(data,round(size(data)./res_mat(:,ll+1)'));
            h5create([save_data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],...
                size(data_temp),'Datatype','uint16','ChunkSize',sub_mat(:,ll+1)');
            h5write([save_data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],data_temp);
        end
    end
end


%% write xml file
docNode = com.mathworks.xml.XMLUtils.createDocument('SpimData');
root = docNode.getDocumentElement;
root.setAttribute('version','0.2');
basePath = docNode.createElement('BasePath');
basePath.setAttribute('type','relative');
basePath.appendChild(docNode.createTextNode('.'));
root.appendChild(basePath);

% sequence description
seqDes = docNode.createElement('SequenceDescription');

imLoader = docNode.createElement('ImageLoader');
imLoader.setAttribute('format','bdv.hdf5');
hdf5 = docNode.createElement('hdf5');
hdf5.setAttribute('type','relative');
hdf5.appendChild(docNode.createTextNode('embryo_data_h5.h5'));
imLoader.appendChild(hdf5);
seqDes.appendChild(imLoader);

viewSets = docNode.createElement('ViewSetups');
Attributes = docNode.createElement('Attributes');
Attributes.setAttribute('name', 'channel');
for ii = 0:view_num-1
    viewSet = docNode.createElement('ViewSetup');
        id = docNode.createElement('id');
        id.appendChild(docNode.createTextNode(num2str(ii)));   % id name
        view_size = docNode.createElement('Size');
        view_size.appendChild(docNode.createTextNode(...
            [num2str(x) ' ' num2str(y) ' ' num2str(z)]));
        voxelSize = docNode.createElement('voxelSize');
        unit = docNode.createElement('unit');
        unit.appendChild(docNode.createTextNode('pixel')); 
        unit_size = docNode.createElement('size');
        unit_size.appendChild(docNode.createTextNode('1 1 1')); 
        voxelSize.appendChild(unit);
        voxelSize.appendChild(unit_size);
        attributes = docNode.createElement('attributes');
            channel = docNode.createElement('channel');
            channel.appendChild(docNode.createTextNode(num2str(ii))); 
        attributes.appendChild(channel); 
    viewSet.appendChild(id);
    viewSet.appendChild(view_size);
    viewSet.appendChild(voxelSize);
    viewSet.appendChild(attributes); 
    viewSets.appendChild(viewSet);

    Channel = docNode.createElement('Channel');
        id = docNode.createElement('id');
        id.appendChild(docNode.createTextNode(num2str(ii)));
        name = docNode.createElement('name');
        if ii == 0
            name.appendChild(docNode.createTextNode('561'));
        else
            name.appendChild(docNode.createTextNode('488'));
        end
    Channel.appendChild(id);
    Channel.appendChild(name);
    Attributes.appendChild(Channel);
end
viewSets.appendChild(Attributes);


seqDes.appendChild(viewSets);

timePoints = docNode.createElement('Timepoints');
timePoints.setAttribute('type','range');
first = docNode.createElement('first');
first.appendChild(docNode.createTextNode('0'));
timePoints.appendChild(first);
last = docNode.createElement('last');
last.appendChild(docNode.createTextNode(num2str(time_num-1)));
timePoints.appendChild(last);
seqDes.appendChild(timePoints);

root.appendChild(seqDes);

% view registration
viewRegs = docNode.createElement('ViewRegistrations');
for tt = 0:time_num-1
    for ii = 0:view_num-1 
        viewReg = docNode.createElement('ViewRegistration');
        viewReg.setAttribute('timepoint',num2str(tt));
        viewReg.setAttribute('setup',num2str(ii));
        viewTrans = docNode.createElement('ViewTransform');
        viewTrans.setAttribute('type','affine');
        affine = docNode.createElement('affine');
        affine.appendChild(docNode.createTextNode(...
            '1 0 0 0 0 1 0 0 0 0 1 0'));
        viewTrans.appendChild(affine);
        viewReg.appendChild(viewTrans);
        viewRegs.appendChild(viewReg);
    end
end
root.appendChild(viewRegs);
xmlwrite([save_data_name '.xml'],docNode);


