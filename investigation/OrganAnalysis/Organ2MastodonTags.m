clc;close all;
clearvars -except movieInfoAll

cell_num = length(movieInfoAll.xCoord);
radius = 11:0.5:14.5;
category = zeros(cell_num,1);
y = 1792;
x = 1818;
z = 253;

% load organ location
path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
name = {'Right Eye', 'Left Eye', 'Brain', 'Somite AB', 'Somite BA', 'Somite AA', 'Somite BB', 'Tail bud'};
for cc = 1:8
    mask = loadchannel(cc);  
    mask = mask > 0;

    cell_im = zeros(y,x,z);
    for ii = sum(movieInfoAll.n_perframe(1:1000))+1 : sum(movieInfoAll.n_perframe)
        cell_loc_y = round((movieInfoAll.yCoord(ii)-1)*2);
        cell_loc_x = round((movieInfoAll.xCoord(ii)-1)*2);
        cell_loc_z = round(movieInfoAll.zCoord(ii)-1);                
        if mask(cell_loc_y, cell_loc_x, cell_loc_z)
            cell_im(cell_loc_y, cell_loc_x, cell_loc_z) = 1;
            category(ii) = cc;
        end
    end
end

function im = loadchannel(channel)
    path = '/work/Mengfan/Embryo/Amat2014/TM1000_reconstruction/export';
    z = 253;
    im = zeros(1792,1818,z);
    for zz = 1:z
        im(:,:,zz) = tifread(fullfile(path, ['SPC0_TM1000_CM0_CM1_CHN00_CHN01.fusedStac-wAnnotations_C' ...
            sprintf('%01d',channel) '_Z' sprintf('%03d',zz-1) '.tif']));
    end
end