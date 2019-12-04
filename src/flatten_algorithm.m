function [ startVolume,startVolumeNormalized,mappedVol] = flatten_algorithm( binaryImg, lambda, rho, relaxFetal, meshParams)
%Flattening Parameterization algorithm.
%Takes as input the placenta segmentation and returns the flattened mesh.
%Outputs the starting mesh generated from the segmentation, and the mesh in
%the flattened space. Both are needed to map MRI intensity values to the
%flattened template.
%Inputs:
%       binaryImage: binary matrix containing the segmentation of the placenta
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
% Algorithm
%Estimate the placenta thickness
surfaceVoxels = binaryImg - imerode(binaryImg, true(3));
distTX = bwdist(surfaceVoxels);
distTX(find(binaryImg == 0)) = -1;
distTX_use = distTX(distTX > 0);
hist_95 = prctile(distTX_use, 95);
rz = hist_95;
%Smooth the image
try
    binaryImg = imgaussfilt3(binaryImg,gaussSigma);
catch
    binaryImg = imgaussfilt3(uint8(binaryImg),gaussSigma);
end
%Create a mesh
[node,elem] = create_mesh(binaryImg,meshParams(1),meshParams(2));
V = elem(:,1:4);
X = node;
V = preprocess_flipVolume(V, X);
%Call the flattening algorithm
[startVolume,startVolumeNormalized, mappedVol] = grad_descent_algo(X,V,lambda,rho,rz, relaxFetal);
end


