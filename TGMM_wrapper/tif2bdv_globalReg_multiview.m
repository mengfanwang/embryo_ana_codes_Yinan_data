function tif2bdv_globalReg_multiview(view_name, tif_folder_path, save_data_name, ...
    timepts_to_process, st_loc, sz_crop, trans_mat)
% view_name: {v} name of views
% tif_folder_path: {v} list of paths
% trans_mat: {t}{v} transformation matrix


%%% scaling setting %%%
if ~isempty(st_loc)
    end_loc = st_loc + sz_crop - 1;
end

%%% parameter %%%
view_num = length(view_name);
res_mat = [1 1 1; 2 2 2;]';
sub_mat = [16 16 16; 8 8 8]';

% write hdf5 data
l_num = size(res_mat,2);

% read tif file
time_num = length(timepts_to_process);
dz_cell = zeros(view_num, 3);
for vv = 1:view_num
    s_ind = char(sprintf("%02d", str2num(view_name{vv})));

%     h5create([save_data_name '.h5'],['/s' s_ind '/resolutions'],size(res_mat));
%     h5write([save_data_name '.h5'], ['/s' s_ind '/resolutions'], res_mat);
%     h5create([save_data_name '.h5'],['/s' s_ind '/subdivisions'],size(sub_mat));
%     h5write([save_data_name '.h5'], ['/s' s_ind '/subdivisions'], sub_mat);
    
    for tt = 1:1 %time_num
        fprintf('Writing time %d\n', tt);
       
        data = tifread(fullfile(tif_folder_path{vv}, timepts_to_process(tt)+'.tif'));
        if ~isempty(st_loc)
            data = data(st_loc(1):end_loc(1), ...
                   st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        end
    
        [x,y,z] = size(data);
        dz_cell(vv, :) = [y x z];
        data = single(data);
    
%         % write h5 file
%         t_ind = char(timepts_to_process(tt));       
%         for ll = 0:l_num-1
%             l_ind = num2str(ll);
%             data_temp = imresize3(data,round(size(data)./res_mat(:,ll+1)'));
%             h5create([save_data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],...
%                 size(data_temp),'Datatype','uint16','ChunkSize',sub_mat(:,ll+1)');
%             h5write([save_data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],data_temp);
%         end
    end
end

% write xml file
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
for vv = 1:view_num
    viewSet = docNode.createElement('ViewSetup');
    id = docNode.createElement('id');
    id.appendChild(docNode.createTextNode(view_name{vv}));   % id name
    view_size = docNode.createElement('Size');
    view_size.appendChild(docNode.createTextNode(...
        [num2str(dz_cell(vv,1)) ' ' num2str(dz_cell(vv,2)) ' ' num2str(dz_cell(vv,3))]));  % not true
    voxelSize = docNode.createElement('voxelSize');
    unit = docNode.createElement('unit');
    unit.appendChild(docNode.createTextNode('pixel')); 
    unit_size = docNode.createElement('size');
    unit_size.appendChild(docNode.createTextNode('1 1 1')); 
    voxelSize.appendChild(unit);
    voxelSize.appendChild(unit_size);
    attributes = docNode.createElement('attributes');
    att_angle = docNode.createElement('angle');
    att_angle.appendChild(docNode.createTextNode(num2str(vv-1)));
    attributes.appendChild(att_angle);
    viewSet.appendChild(id);
    viewSet.appendChild(view_size);
    viewSet.appendChild(voxelSize);
    viewSet.appendChild(attributes);
    viewSets.appendChild(viewSet);
end
Attributes = docNode.createElement('Attributes');
Attributes.setAttribute('name','angle');
for vv = 1:view_num
    angle = docNode.createElement('Angle');
    angle_id = docNode.createElement('id');
    angle_id.appendChild(docNode.createTextNode(num2str(vv-1)));
    angle_name = docNode.createElement('name');
    angle_name.appendChild(docNode.createTextNode(num2str(90*(vv-1))));
    angle_axis = docNode.createElement('axis');
    angle_axis.appendChild(docNode.createTextNode('0.0 1.0 0.0'));
    angle_degreees = docNode.createElement('degrees');
    angle_degreees.appendChild(docNode.createTextNode(num2str(90*(vv-1))));

    angle.appendChild(angle_id);
    angle.appendChild(angle_name);
    angle.appendChild(angle_axis);
    angle.appendChild(angle_degreees);
    Attributes.appendChild(angle);
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

% % view registration
% trans_mat = cell(1,time_num);
% for tt = time_num:-1:1   % registred from tt to tt+1
%     if tt == time_num
%         trans_mat{time_num} = eye(4,4);
%         trans_mat{time_num}(3,3) = 5.412; % temp for Joaquin data
%     else
%         trans_mat{tt} = tform{tt}.T*trans_mat{tt+1};
%     end
% end
viewRegs = docNode.createElement('ViewRegistrations');
for tt = 0:time_num-1
    for vv = 1:view_num 
        trans_str = cell(1,12);
        trans_tmp = trans_mat{tt+1}{vv}';
        for jj = 1:12
            trans_str{jj} = num2str(trans_tmp(jj));
        end

        viewReg = docNode.createElement('ViewRegistration');
        viewReg.setAttribute('timepoint',num2str(tt));
        viewReg.setAttribute('setup',view_name{vv});
        viewTrans = docNode.createElement('ViewTransform');
        viewTrans.setAttribute('type','affine');
        affine = docNode.createElement('affine');
%         affine.appendChild(docNode.createTextNode(...
%             '1 0 0 0 0 1 0 0 0 0 1 0'));
        affine.appendChild(docNode.createTextNode(...
            join(trans_str, ' ')));
        viewTrans.appendChild(affine);
        viewReg.appendChild(viewTrans);
        viewRegs.appendChild(viewReg);
    end
end
root.appendChild(viewRegs);
xmlwrite([save_data_name '.xml'],docNode);


