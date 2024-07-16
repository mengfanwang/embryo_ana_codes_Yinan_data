clc;clear;close all;

path = '/work/Mengfan/paper_fig/Amat2014_visual';
numFrames = 50;
video = zeros(740,740*2,3,numFrames);

file1 = 'CartesianMapping_769.avi';
vidObj = VideoReader(fullfile(path, file1));

for ii = 1:numFrames
    for jj = 1:10
        vidFrame = readFrame(vidObj);
    end
    video(:,1:740,:,ii) = vidFrame(100:839,135:874,:);
end

file2 = 'TGMM_720.avi';
vidObj = VideoReader(fullfile(path, file2));
for ii = 1:numFrames
    for jj = 1:10
        vidFrame = readFrame(vidObj);
    end
    video(:,741:end,:,ii) = vidFrame(100:839,135:874,:);
end


%% write as gif
filename = fullfile(path,'compare_500.gif');
% write the frames to the video
for i=1:numFrames
    frame = video(:,:,:,i)/255;  
    frame = insertText(frame,[660 35], sprintf('T = %0.1fh', i/10+4), 'FontSize', 40, 'BoxColor', [1 1 1]);
    [A,map] = rgb2ind(frame,256);
    if i == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.2);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.2);
    end
end