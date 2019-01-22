function [ binaryImg_cleaned ] = preprocess_seg_holes( binaryImg)
%This function loads in a segmentation image file (3D matrix) and fills
%holes and removes floating components.
%Inputs:
%   binaryImg: 3D binary matrix. 
%Outputs:
%   binaryImg_cleaned: 3D binary matrix with the largest segmentation
%   component, and holes cleared.

CC = bwconncomp(binaryImg,26);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
binaryImg(:) = 0;
binaryImg(CC.PixelIdxList{idx}) = 1;
binaryImg_cleaned = imfill(binaryImg,'holes');
end

