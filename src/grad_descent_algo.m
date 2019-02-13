function [ startVolume, startVolumeNormalized, mappedVolume] = grad_descent_algo( X, T, lambda, rho, rz )
%Mapping algorithm optimized using gradient descent. Parameterizes the
%placenta by mapping to a flattened space
%Inputs: 
%       X: mesh vertices
%       T: mesh triangulation
%       lambda: regularization parameter
%       rho: line search descent factor
%       rz: template half-height
%       savePath: location to save output files
%Outputs:
%       startVolume: volumetric mesh of the placenta in the original image
%           space
%       startVolumeNormalized: volumetric mesh of the placenta in the
%       original image space, centered at 0 and rotated by moments of
%       inertia
%       mappedVolume: volumetric mesh of the placenta in the flattened space. 
close all;
%optimization parameters
flattenMaternalOnly = 0;
flattenMedialAxis = 0;
medialAxisHeight = 1;
numIterations =75001; 
gradThresh = 1e-4;
objThresh = 400;
normThresh = 150;
lineSearchLim = 200;
useRim = 1;

%normalize and PCA
if(size(X,2)>size(X,1))
    X = X';
end
Xorig = X;
meanX = mean(X);
X = X - meanX;
[coeff, scores, latent]  = pca(X);
X = X*coeff;

nPts = length(T);
DT=triangulation(T,X);
T = DT.ConnectivityList;
X0P = DT.Points;
X0 = gpu_generate_tets(T, X0P);
Xvol = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
Vol = sum(Xvol);

D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];
O = pagefun(@mtimes, X0,D);
Oinv = pagefun(@inv, O);
OT = pagefun(@transpose, O);
OinvT = pagefun(@transpose, Oinv);
OTO = pagefun(@mtimes, OT,O);

startVolume = triangulation(T,Xorig);
%extract the boundary triangles
FBtri = freeBoundary(DT);
boundaryVertices = unique(FBtri);

%compute the voronoi areas per boundary node.
[Ts, Xs] = freeBoundary(DT);
voronoiAreas = voronoi_area(Ts, unique(Ts), Xs);
voronoiAreas = voronoiAreas / sum(voronoiAreas);
%normalize volume
Xvol = Xvol /Vol;

%Label maternal and fetal sides, and find the placenta rim.
[Ts, Xs] = freeBoundary(DT);
[binRim, binNorth, binSouth, binNorthNoRim, binSouthNoRim] = extract_rim_spectral_clustering(Ts, Xs);
binMap = map_boundary_node_index(FBtri, Ts);
dist_1 = gpu_compute_distortion_rim(X0P, voronoiAreas, rz, binNorth, binSouth, binMap);
dist_2 = gpu_compute_distortion_rim(X0P, voronoiAreas, rz, binSouth, binNorth, binMap);
%flip, to correctly chick north and south hemisphere allocation.
if(dist_2 <= dist_1)
    binNT = binNorth;binST = binSouth;
    binNorth = binST;binSouth = binNT;
    binNT = binNorthNoRim; binST = binSouthNoRim;
    binNorthNoRim = binST; binSouthNoRim = binNT;
end
[maternalNodes, fetalNodes] = label_maternal_fetal(Ts, Xs, binNorth, binSouth);
%flip, so fetal side is at the top.
if(isequal(maternalNodes, binNorth))
    Rflip = [1 0 0;0 1 0;0 0 -1];
    X0P = X0P*Rflip;
    binNT = binNorth;
    binST = binSouth;
    binNorth = binST;
    binSouth = binNT;
end

%Normalized start volume
startVolumeNormalized = triangulation(T,X0P);

%Gradient Preconditioning (not used)
L = -cotmatrix(X0P,T);
nv = length(X0P);
%add the null space back
nullSpaceBasis = ones(1,nv)./sqrt(nv);
L = L+nullSpaceBasis'*nullSpaceBasis;
AcotL = chol(L);

%% Gradient Descent
gradNorms = zeros(1,numIterations);
countFunc = 0;
countNormFunc = 0;
ii=1;
brokeLineSearch = 0;
converged = 0;

%Start timing
tic
%initialization
XPTemp = X0P;
XP = X0P;
isStuck = 0;
brokeLineSearch = 0;


X = gpu_generate_tets(T, XPTemp);
%Compute the gradient
grad = gpu_compute_grad(X,T,lambda, Xvol, Oinv, OinvT, OTO);
grad = gpu_compute_grad_rim(grad, XP, voronoiAreas, rz, binNorth, binSouth, binMap);
gradNorms(ii) =gather(norm(grad,'fro'));
%Gradient descent
while(norm(grad,'fro')>gradThresh && ii<numIterations && countFunc<objThresh && isStuck == 0 && brokeLineSearch == 0)
    count = 0;
    X = gpu_generate_tets(T, XPTemp);
    distData = gpu_compute_distortion_rim(XP, voronoiAreas, rz, binNorth, binSouth, binMap);
    origDist = gpu_compute_distortion(X, lambda, Xvol, Oinv);
    origDist = origDist + distData;
    
    grad = gpu_compute_grad(X,T,lambda, Xvol, Oinv, OinvT, OTO);
    grad = gpu_compute_grad_rim(grad, XP, voronoiAreas, rz, binNorth, binSouth, binMap);
    if(ii == 1)
        gradO = grad;
    end
    gradX = gpu_generate_tets(T,grad);
    
    %Compute step size to ensure injectivity
    eta = real(gpu_compute_minT( X,gradX ));
    eta = eta*0.9;
    orig_eta =eta;
    XPTemp = XP - eta*grad;
    
    %check if actually injective.
    distVolume = compute_distortion_J(T, T, XP, XPTemp);
    if(min(distVolume)<=0)
        while( min(distVolume)<=0)
            eta =  eta/10;
            XPTemp = XP - eta*grad;
            distVolume = compute_distortion_J(T, T, XP, XPTemp);
        end
    end
    
    X = gpu_generate_tets(T, XPTemp);
    [distData, ~, ~] = gpu_compute_distortion_rim(XPTemp, voronoiAreas, rz, binNorth, binSouth, binMap);
    [distJTemp, ~, ~] = gpu_compute_distortion(X, lambda, Xvol, Oinv);
    distI = distJTemp + distData;
    
    gradNorms(ii) = gather(norm(grad,'fro'));
    %Line Search
    while(distI >= origDist)%+beta*eta*norm(grad,'fro'))
        eta = eta*rho;
        XPTemp = XP - eta*grad;
        X = gpu_generate_tets(T, XPTemp);
        [distData, ~,~] = gpu_compute_distortion_rim(XPTemp, voronoiAreas, rz, binNorth, binSouth, binMap);
        [distJTemp,~,~] = gpu_compute_distortion(X, lambda, Xvol, Oinv);
        distI = distJTemp + distData;
        count = count+1;
        if(count>lineSearchLim)
            brokeLineSearch = 1;
            fprintf ('Iteration: %d,  grad: %d, obj: %d \n', ii, gradNorms(ii),distI)
            fprintf('Took too many line search loops at iteration %d \n',ii);
            break;
        end        
    end
    ii = ii+1;
    XP = XPTemp;
    gradNorms(ii) = gather(norm(grad,'fro'));  
    if(ii<100)
        fprintf ('Iteration: %d, grad: %d, obj: %d \n', ii, gradNorms(ii),distI)
        fprintf('\n')
        %save([savePath,'/startVolume'],'startVolume','-v7.3');
        %save([savePath,'/startVolumeNormalized'],'startVolumeNormalized','-v7.3');
    end
    if(mod(ii,100)==0)
        fprintf ('It: %d, grad: %d, obj: %d \n', ii, gradNorms(ii),distI)
        fprintf('\n')
    end
    if(abs(distI-origDist)<1e-6*lambda)
        countFunc = countFunc+1;
    else
        countFunc = 0;
    end
    
    if(eta < 1e-50)
        isStuck = 1;
        fprintf ('It: %d, Orig eta: %d, eta: %d, grad: %d, obj: %d \n', ii, orig_eta, eta, gradNorms(ii),distI)
        break; %algorithm cannot move
    end
end
%% Output data
Xtrans = gather(XP)+ meanX;
XP = Xtrans;
mappedVolume=triangulation( T,gather(XP));
toc

