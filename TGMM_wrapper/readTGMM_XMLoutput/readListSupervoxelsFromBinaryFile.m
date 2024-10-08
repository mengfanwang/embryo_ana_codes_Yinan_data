% Reads binary file written with C++ function lineageHyperTree::writeListSupervoxelsToBinaryFile
% 
% Each XML file output by TGMM has a file with the same name but with the .svb extenstion (from Super-Voxel Binary file). This a custom-created binary file which stores all the indexes in the image for each super-voxel at a specific time point. In this way, and using svIdx we can recover the precise segmentation of each tracked object. 
% In order to do that, you need to to use the Matlab function [svList, sizeIm] = readListSupervoxelsFromBinaryFile( filename )
% 
% INPUT:
% 
% filename:		string containing the path to the svb file generated by TGMM code. For exmaple, if the TGMM oputput is located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht and you want to read time point 100, then filenameXML = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame0100.svb'
% 
% OUTPUT:
% 
% svList:		cell array of length S, where S is the number of super-voxels generated at the time point specified by filename. svList{i} contains the indexes of all the voxels in the image belonging to the i-th supervoxel in ascending order. This list is exactly the same as the variable PixelIdxList returned by Matlab functions such as regionprops or bwconncomp.
% 		Use svList{svIdx+1} in order to access the correct super-voxel, since svIdx follow C-indexing convention.
% 
% sizeIm:		1x3 integer vector containing the size of the image stack along each dimension. This is useful to recover (x,y,z) positions from svList indexes using the function ind2sub.

function [svList, sizeIm] = readListSupervoxelsFromBinaryFile( filename )



fid = fopen(filename, 'rb');

if( fid < 0 )
   svList = [];
   sizeIm = [];
   return;
end

numSv = fread(fid,1,'int32');
svList = cell( numSv, 1);

for kk = 1: numSv
   svStruct = readSupervoxelFromBinary( fid ); 
   svList{kk} = svStruct.PixelIdxList;
end

sizeIm = svStruct.dataDims';
fclose(fid);