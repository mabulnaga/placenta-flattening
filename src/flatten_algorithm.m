function [ startVolume,startVolumeNormalized,mappedVol] = flatten_algorithm( binaryImg, lambda, rho, relaxFetal, meshParams,varargin)
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
%       varargin:
%           img_spacing. By default, code assumes isotropic spacing.
%           Otherwise, code inputs a 3x3 diagonal matrix, with the image
%           spacing to rescale the mesh.
%Outputs:
%       startVolume: mesh generated from the placenta mask
%       startVolumeNormalized: startVolume centered at the origin and
%       rotated by moments of inertia.
%       mappedVol: mapped mesh in the flattened space.

%Smoothing parameter
gaussSigma = 1.5;
isotropic = 1;
if(nargin > 5)
    isotropic = 0;
    img_spacing_mat = varargin{1};
end
% Algorithm
%Estimate the placenta thickness
if(isotropic == 1)
    surfaceVoxels = binaryImg - imerode(binaryImg, true(3));
    distTX = bwdist(surfaceVoxels);
    distTX(find(binaryImg == 0)) = -1;
    distTX_use = distTX(distTX > 0);
    hist_95 = prctile(distTX_use, 95);
    rz = hist_95;
end 
%Smooth the image
try
    binaryImg = imgaussfilt3(binaryImg,gaussSigma);
catch
    binaryImg = imgaussfilt3(uint8(binaryImg),gaussSigma);
end
%Create a mesh
[node,elem] = create_mesh(binaryImg,meshParams(1),meshParams(2));
T = elem(:,1:4);
X = node;
T = preprocess_flipVolume(T, X);
% get the r_z estimated from the mesh.
if(isotropic == 0)
    %first, extract the boundary mesh.
    X =X*img_spacing_mat; 
    Ts = freeBoundary(triangulation(T,X));
    % now, for every interior vertex, find the closest euclidean distance
    % to the boundary.
    interior_vertices =  X;
    interior_vertices(unique(Ts),:) = [];
    vec_distances = zeros(1,length(interior_vertices));
    boundary_vertices = X(unique(Ts),:);
    for i = 1 : length(interior_vertices)
        v = interior_vertices(i,:);
        v = repmat(v,length(boundary_vertices),1);
        dist = sqrt(sum((v - boundary_vertices).^2,2));
        vec_distances(i) = min(dist);
    end
    rz = prctile(vec_distances,95);
end
%Call the flattening algorithm
[startVolume,startVolumeNormalized, mappedVol] = grad_descent_algo(X,T,lambda,rho,rz, relaxFetal);
end


