function [grad, gradNorth, gradSouth] = gpu_compute_grad_rim(grad, P, surfaceScale, rz, northHem, southHem, binMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(isempty(northHem) && isempty(southHem))
   gradNorth = [];
   gradSouth = [];
   return;
end
northPts = P(binMap(northHem),:);
northAreas = surfaceScale(northHem);
southPts = P(binMap(southHem),:);
southAreas = surfaceScale(southHem);

northDist = (northPts(:,3) - rz);
northDist = 2*northDist.*northAreas;
gradNorth = northDist;

southDist = (southPts(:,3) + rz);
southDist = 2*southDist.*southAreas;
gradSouth = southDist;

grad(binMap(northHem),3) = grad(binMap(northHem),3) + gradNorth;
grad(binMap(southHem),3) = grad(binMap(southHem),3) + gradSouth;

end

