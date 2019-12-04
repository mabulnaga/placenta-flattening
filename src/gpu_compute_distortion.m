function [ totalD, distJMax, distJnode] = gpu_compute_distortion( X, lambda, Xvol,Oinv)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
%X0volumes = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
%Xvol = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
J = gpu_compute_J(X,Oinv);
JT = pagefun(@transpose, J);
JJT =pagefun(@mtimes,J,JT);
Jinv = pagefun(@inv, J);
JTinv = pagefun(@inv, JT);
JJTinv = pagefun(@mtimes, Jinv, JTinv);

diags = reshape(bsxfun(@plus,(0:size(J,3)-1)*9,(1:4:9).'),3,1,[]);
j = JJT(diags);
ji = JJTinv(diags);
%D1 = arrayfun(@(x) gather(trace(JJT(:,:,x))), 1:size(JJT,3));
%D2 = arrayfun(@(x) gather(trace(JJTinv(:,:,x))), 1:size(JJTinv,3));
%D = D1 + D2;
D = sum(j)+sum(ji);
D = reshape(D,1,length(D));
scaledD = gpuArray(Xvol).*D;
%alpha is the weight on regularization (can think of it as lambda)
totalD = lambda * sum(scaledD);
DJ = totalD;
distJnode = scaledD;
%scaledD: scaled distortion (by tet volume), per tet (vector).
%totalD: summed distortion
%Convert to GPU...

%sphere norms
% boundaryPts = P(uniqueVerts,:);
% norms = sum(boundaryPts.^2,2);
% norms = (norms - r^2).^2;
% totalD = totalD + alpha*sum(norms);
% DSphere = alpha*sum(norms);

%ellipsoid norms
% boundaryPts = P(uniqueVerts,:);
% boundaryPts(:,1) = boundaryPts(:,1)/r;
% boundaryPts(:,2) = boundaryPts(:,2)/ry;
% boundaryPts(:,3) = boundaryPts(:,3)/rz;
% norms = sum(boundaryPts.^2,2);
% norms = (norms - 1).^2;
% DSphere = sum(surfaceScale.*norms);


%for statistics, also output the worst case ones for the data distortion
%and regularization distortiton
[distJMax, distJMaxInd] = max(scaledD);
distJMax = [distJMax; distJMaxInd];
end

