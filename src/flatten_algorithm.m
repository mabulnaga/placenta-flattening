function [ startVolume,startVolumeNormalized,mappedVol] = flatten_algorithm( segmentationMask, lambda, rho, relaxFetal, meshParams)
%Flattening Parameterization algorithm.
%Takes as input the placenta segmentation and returns the flattened mesh.
%Outputs the starting mesh generated from the segmentation, and the mesh in
%the flattened space. Both are needed to map MRI intensity values to the
%flattened template.
%Inputs:
%       segmentationMask: binary matrix containing the segmentation of the placenta
%       lambda: regularization parameter governing the tradeoff between
%       template fit and shape distortion
%       rho: line search descent parameter
%       relaxFetal: relax fetal side: flag to indicate whether to relax the
%               fetal side of the placenta after optimization (1) or not
%               (0). Default: 0.
%       meshParams: vector with 2 elements that contains the meshing
%       parameters. See iso2mesh documentation for additional details.
%Outputs:
%       startVolume: mesh generated from the placenta mask
%       startVolumeNormalized: startVolume centered at the origin and
%       rotated by moments of inertia.
%       mappedVol: mapped mesh in the flattened space.


%Smoothing parameter
gaussSigma = 1.5;
% dilation parameter
dilateNum = 3;
% smooth and dilate the mask
segmentationMask = smooth_segmentation(segmentationMask, dilateNum, gaussSigma);
%Estimate the placenta thickness
rz = estimate_placenta_thickness(segmentationMask);

%Create a mesh
[node,elem] = create_mesh(segmentationMask,meshParams(1),meshParams(2));
V = elem(:,1:4);
X = node;
V = preprocess_flipVolume(V, X);
%Call the flattening algorithm
[startVolume,startVolumeNormalized, mappedVol] = grad_descent_algo(X,V,lambda,rho,rz, relaxFetal);
end


