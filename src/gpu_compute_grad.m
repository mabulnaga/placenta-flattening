function [ grad ] = gpu_compute_grad( X,T,lambda, Xvol,Oinv, OinvT, OTO)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];
%X0volumes = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
%Xvol = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));   
%global O Oinv OT OinvT OTO
%O = pagefun(@mtimes, X0,D);
E = pagefun(@mtimes,X,D);
%Oinv = pagefun(@inv, O);
%OT = pagefun(@transpose, O);
XD = E;
XT = pagefun(@transpose, X);
%OinvT = pagefun(@transpose, Oinv);
% EinvT = pagefun(@inv, E);
% EinvT = pagefun(@transpose, EinvT);
% [scaledD, totalD] =gpu_compute_distortion(X, X0);

 t1 = pagefun(@mtimes, XD, Oinv);
 t1 = pagefun(@mtimes, t1, OinvT);
 t1 = pagefun(@mtimes, t1, D');
 t1 = 2*t1;
% 
  t22 = pagefun(@mtimes, D',XT);
  t22 = pagefun(@mtimes, t22, X);
  t22 = pagefun(@mtimes, t22, D);
  t22 = pagefun(@inv, t22);
  
  %t23 = pagefun(@mtimes, OT,O);
  t23 = OTO;
  t2 = 2*XD;
  t2 = pagefun(@mtimes, t2, t22);
  t2 = pagefun(@mtimes, t2, t23);
  t2 = pagefun(@mtimes, t2, t22);
  t2 = pagefun(@mtimes, t2, D');
% t22 = pagefun(@mtimes, OinvT, D');
% t22 = pagefun(@mtimes, t22, XT);
% t22 = pagefun(@mtimes, t22, X);
% t22 = pagefun(@mtimes, t22, D);
% t22 = pagefun(@mtimes, t22, Oinv);
% t22 = pagefun(@inv, t22);
% t22 = pagefun(@mtimes, t22, t22);
% 
% t23 = pagefun(@mtimes, OinvT, D');
% 
% t2 = pagefun(@mtimes, t21, t22);
% t2 = pagefun(@mtimes, t2, t23);

grad = t1 - t2;
grad = bsxfun(@times, grad, reshape(Xvol,[1,1,size(grad,3)]));
%multiply by the regularization term penalty
%everything above here was the first term of the summation, now will add
%the derivative of the volume term.
%Now, convert to per vertex gradient.
T = reshape(T',1,[]);
gradMat = reshape(grad,3,[]);
z1 = accumarray(T',gradMat(1,:)')';
z2 = accumarray(T',gradMat(2,:)')';
z3 = accumarray(T',gradMat(3,:)')';
gradientA = lambda* [z1;z2;z3]';

%Now, need gradient for the spherical constraint term.
% boundaryPts = P(uniqueVerts,:);
% norms = sum(boundaryPts.^2,2);
% norms = 4*(norms-r^2);
% norms = [norms , norms, norms];
% gradSphere = alpha*norms.*boundaryPts;
%Ellipsoid gradient
% boundaryPts = P(uniqueVerts,:);
% boundaryPts(:,1) = boundaryPts(:,1)/r;
% boundaryPts(:,2) = boundaryPts(:,2)/ry;
% boundaryPts(:,3) = boundaryPts(:,3)/rz;
% norms = sum(boundaryPts.^2,2);
% norms = 4*(norms - 1);
% norms = surfaceScale.*norms;
% norms = [norms/r, norms/ry, norms/rz];
% gradSphere = norms.*boundaryPts;

%gradientA(uniqueVerts,:) = gradientA(uniqueVerts,:)+gradSphere;
grad = gradientA;
end

