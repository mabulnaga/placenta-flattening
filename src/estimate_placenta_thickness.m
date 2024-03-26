function [rz] = estimate_placenta_thickness(segmentationMask)
%Estimates the thickness of the placenta
% Arguments:
%           segmentationMask: 3D binary image with placenta segmentation
% Returns:
%           rz: (scalar) estimated z-axis thickness.

% estimate the voxels on the surface
surfaceVoxels = segmentationMask - imerode(segmentationMask, true(3));
distTX = bwdist(surfaceVoxels);
distTX(find(segmentationMask == 0)) = -1;
distTX_use = distTX(distTX > 0);

% this is another method to get the thickness, but it
% overestimates. Only call when the previous one gives nans.
% May also want to do more experiments with this one, though.
if(isempty(distTX_use))
    distTX = bwdist(~logical(segmentationMask));
    distTX(find(segmentationMask<0.05))=-1;
    distTX_use = distTX(distTX>0);
end

hist_95 = prctile(distTX_use, 95);
rz = hist_95;

end

