function [ detJ ] = compute_distortion_J(T0, T, X0P, XP )
%computes the determinant of the Jacobian matrix of a transformation from
%original set of points X0P to a new set, XP. Returns the determinant of
%the Jacobian per tet.
%Inputs: T0 original connectivity list
%       T: new connectivity list (should be same as T0)
%       X0P: original set of poitns
%       XP: new set of pts
%Output: detJ, size 1xnT, determinant of jacobian per tet


X0 = gpu_generate_tets(T0,X0P);
X = gpu_generate_tets(T, XP);
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
O = pagefun(@mtimes, X0,D);
Oinv = pagefun(@inv, O);

J = gpu_compute_J(X, Oinv);

A = J;
a = A(1,1,:);b =A(1,2,:);c = A(1,3,:);d = A(2,1,:);e = A(2,2,:); 
f= A(2,3,:); g = A(3,1,:); h = A(3,2,:); i= A(3,3,:);
  
detJ= a.*(e.*i - f.*h) - b.*(d.*i - f.*g) + c.*(d.*h - e.*g);
end
% XnewVol = arrayfun(@(x) tet_volume(X(:,:,x)), 1:size(X,3));
% XoVol =  arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
% x = X(:,:,2);
% z = X0(:,:,2);
% Jj = (x*D)*inv(z*D);

