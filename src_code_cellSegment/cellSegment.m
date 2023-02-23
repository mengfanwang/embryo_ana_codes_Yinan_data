function [refine_res, eig2d, eig3d, threshold_res] = cellSegment(org_im, repeat_flag)
if nargin == 1
    repeat_flag = true;
end
if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
end

addpath('src_code_matlab');
addpath('src_code_cellSegment');
addpath('src_code_cellTracker');
addpath('./');


minIntensity = 0;
%% synQuant for the 40 frames
% add synQuant java path
Pij = fullfile('src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
%p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.4.jar');
p0 = fullfile('src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);%

q.minIntensity = minIntensity;
sigma = 1;
sm_im = imgaussfilt3(org_im,sigma);
%q.posEigMap = eig_res_3d{i}>0;
% 3D version
[zMap, synId, fMap] = Synquant4Embryo_Paramater(sm_im, q);
synId = uint16(synId);
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);

%% refine results from synQuant

sigma = 3; % 2d principal curvature smooth scale is smaller
fMaps = org_im > minIntensity;
[eig2d, ~] = principalCv2d(org_im, synId, sigma, fMaps);
sigma = [5,1];
[eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps);

%% calculate the variance map of all frames
scale_term = 5000;

vid = 255*org_im/scale_term;
varMap = cell(3,2);
[varMap{1,1}, varMap{2,1},varMap{3,1}] = ...
    calVarianceStablizationBY(vid, 0.8, 3);
vid_stb = sqrt(vid+3/8);
[varMap{1,2}, varMap{2,2},varMap{3,2}] = ...
    calVarianceStablizationBY(vid_stb, 0.8, 3);

%% region refine based on 4d information (infor across >1 frame)
in_im = vid;
eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
eigAll{1} = eig2d;
eigAll{2} = eig3d;
varMapAll = varMap;

[newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
    eigAll,in_im, varMapAll, []);%, i

refine_res = uint32(newIdMap);
threshold_res = uint8(thresholdMap);
if repeat_flag == false
    return;
end
%% second time of synQuant
% add synQuant java path
Pij = fullfile('../../ParticleTracking', 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('../../ParticleTracking','src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
%p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.4.jar');
p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);%

sigma = 1;
sm_im = imgaussfilt3(org_im,sigma);
posEigMap = eig3d>0;
% 3D version
id_mat_2nd = Synquant4Embryo_2iter(sm_im, refine_res, posEigMap);

%% refine the new detected cells
id_mat_com = refine_res;

cell_append = max(id_mat_2nd(:));
valid_fg = id_mat_com>0;
id_mat_com(valid_fg) = id_mat_com(valid_fg) + cell_append;
id_mat_com = id_mat_com + id_mat_2nd;

refine_res_1st = refine_res;
threshold_res_1st = threshold_res;
%svf = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\segmentation_from_synQuant\';

synId = id_mat_com;

test_ids = 1:double(max(id_mat_2nd(:)));
[newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
    eigAll,vid, varMapAll, test_ids, []);

valid_fg = refine_res_1st>0;
refine_res_1st(valid_fg) = refine_res_1st(valid_fg) + ...
    uint32(max(newIdMap(:)));
refine_res = uint32(newIdMap) + refine_res_1st;
threshold_res = uint8(thresholdMap) + threshold_res_1st;

