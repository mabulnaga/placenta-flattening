function [ startVolume, startVolumeNormalized, mappedVolume] = grad_descent_algo( X, T, lambda, rho, rz, relaxFetal )
%Mapping algorithm optimized using gradient descent. Parameterizes the
%placenta by mapping to a flattened space
%Inputs:
%       X: mesh vertices
%       T: mesh triangulation
%       lambda: regularization parameter
%       rho: line search descent factor
%       rz: template half-height
%       relaxFetal: flag to relax the fetal side after convergence
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

%parameters to govern when to start optimizing the template height
rzStart = gradThresh/10;
rzEnd = gradThresh;
optimizationCount = 2;

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
gradNorms = zeros(1,numIterations,'gpuArray');
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
ignoreNodes = [];

X = gpu_generate_tets(T, XPTemp);
%Compute the gradient
grad = gpu_compute_grad(X,T,lambda, Xvol, Oinv, OinvT, OTO);
grad = gpu_compute_grad_rim(grad, XP, voronoiAreas, rz, binNorth, binSouth, binMap);
gradNorm = norm(grad,'fro');
gradNorms(ii) = gradNorm; %prev. gather
%Gradient descent
for optimizationCounter = 1 : optimizationCount
    while(gradNorm>gradThresh && ii<numIterations && countFunc<objThresh && isStuck == 0 && brokeLineSearch == 0)
        count = 0;
        X = gpu_generate_tets(T, XPTemp);
        distData = gpu_compute_distortion_rim(XP, voronoiAreas, rz, binNorth, binSouth, binMap);
        origDist = gpu_compute_distortion(X, lambda, Xvol, Oinv);
        origDist = origDist + distData;
        
        %compute the gradients
        grad = gpu_compute_grad(X,T,lambda, Xvol, Oinv, OinvT, OTO);
        [grad,~,~,grad_rz] = gpu_compute_grad_rim(grad, XP, voronoiAreas, rz, binNorth, binSouth, binMap);
        grad(ignoreNodes,:) = 0;
        gradNorm = norm(grad,'fro');
        grad = grad/gradNorm;
        gradX = gpu_generate_tets(T,grad);
        
        %Compute step size to ensure injectivity
        eta = real(gpu_compute_minT( X,gradX ));
        eta = eta*0.9;
        orig_eta =eta;
        
        %set the step size for the height optimization
        if(optimizationCounter > 1)
            eta = min(eta,10);
            eta_rz = eta;
        else
           eta_rz = 0; 
        end
        
        %move by the gradient
        XPTemp = XP - eta*grad;
        rz_temp = rz - eta_rz*grad_rz;

        
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
        [distData, ~, ~] = gpu_compute_distortion_rim(XPTemp, voronoiAreas, rz_temp, binNorth, binSouth, binMap);
        [distJTemp, ~, ~] = gpu_compute_distortion(X, lambda, Xvol, Oinv);
        distI = distJTemp + distData;
        
        gradNorms(ii) = gradNorm;
        %Line Search
        while(distI >= origDist)%+beta*eta*norm(grad,'fro'))
            eta = eta*rho;
            eta_rz = eta_rz*rho;
            XPTemp = XP - eta*grad;
            rz_temp = rz - eta_rz * grad_rz;
            X = gpu_generate_tets(T, XPTemp);
            [distData, ~,~] = gpu_compute_distortion_rim(XPTemp, voronoiAreas, rz_temp, binNorth, binSouth, binMap);
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
        gradNorms(ii) = gradNorm;
        
        %store the new vertex, rz locations.
        XP = XPTemp;
        rz = rz_temp;
        
        if(ii<100)
            fprintf ('Iteration: %d, grad: %d, obj: %d \n', ii, gradNorms(ii),distI)
            fprintf('\n')
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
end
%% Output data
Xtrans = gather(XP)+ meanX;
XP = Xtrans;
mappedVolume=triangulation( T,gather(XP));
toc

