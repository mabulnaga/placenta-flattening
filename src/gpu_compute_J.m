function [ J ] = gpu_compute_J( X,Oinv)
%gpu_compute_J Summary of this function goes here
%   Assumes X, X0 are GPU arrays. Computes the Jacobian as a 3D-array,
%   where the third dimension is the Jacobian of the tetrahedron.
%   Convention: X is 3x4.
D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];
%O = pagefun(@mtimes, X0,D);
E = pagefun(@mtimes,X,D);
%inverses = pagefun(@inv, O);
J = pagefun(@mtimes, E,Oinv);

end

