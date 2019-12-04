function [grad, gradNorth, gradSouth, grad_rz] = gpu_compute_grad_rim(grad, P, surfaceScale, rz, northHem, southHem, binMap)
%computes the gradient on the template term when mapping to two planes.
%Also outputs the gradient for the height template.
%Also, outputs the optimal r_z value.
if(isempty(northHem) && isempty(southHem))
   gradNorth = [];
   gradSouth = [];
   grad_rz = 0;
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

grad_rz = -1*sum(northDist) + sum(southDist);
end