function [rimDist, rimDistMax, rimDistDistribution] = gpu_compute_distortion_rim(P, surfaceScale, rz, northHem,southHem, binMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%rimDistDistribution: -1 on the rim, otherwise positive voxel-wise distance
%to template. always positive, and NOT squared distance. also, NOT weighted
%by the surface Area.
if(isempty(northHem) && isempty(southHem))
   rimDist = 0; 
   rimDistMax =[0;0];
   rimDistDistribution = zeros(size(binMap));
   return
end
rimDistDistribution = NaN*ones(length(surfaceScale),1);
northPts = P(binMap(northHem),:);
northAreas = surfaceScale(northHem);
southPts = P(binMap(southHem),:);
southAreas = surfaceScale(southHem);

northDist = (northPts(:,3) - rz).^2;
rimDistDistribution(northHem) = gather(northDist.^(1/2)); %need to be positive
northDist = northDist.*northAreas;
northDistSum = sum(northDist);

southDist = (southPts(:,3) + rz).^2;
rimDistDistribution(southHem) = gather(southDist.^(1/2));
southDist = southDist.*southAreas;
southDistSum = sum(southDist);

rimDist = southDistSum + northDistSum;

allHems = [northHem, southHem];
[maxDist, maxDistInd] = max([northDist;southDist]);
maxDist = gather(maxDist);
maxDistInd = binMap(allHems(maxDistInd));

rimDistMax = [maxDist;maxDistInd];
if(isempty(rimDistMax))
    rimDistMax=[0;0];
end
end

