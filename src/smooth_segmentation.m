function [segmentationMask] = smooth_segmentation(segmentationMask, dilateNum, gaussSigma)
%Function to dilate and smooth a segmentation mask.
%Order of operations: dilate then smooth.
%Args:
%     segmentationMask: 3D binary image with placenta segmentation
%     dilateNum: (scalar) parameter indicating how much to dilate mask by.
%     gaussSigma: (scalar) gauss filtering smoothing variance
%Returns:
%     segmentationMask: smoothed and dilated mask.

if(dilateNum ~=0)
    segmentationMask = imdilate(segmentationMask,true(dilateNum));
    segmentationMask(segmentationMask>0.1) = 1;
end            
%Smooth the image
try
    segmentationMask = imgaussfilt3(segmentationMask,gaussSigma);
catch
    segmentationMask = imgaussfilt3(uint8(segmentationMask),gaussSigma);
end

end

