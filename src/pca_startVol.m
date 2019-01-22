function [X] = pca_startVol(X)
%Normalizes the start volume points and rotates by moments of inertia.
%Input: vector of points, size nvx3. Output: rotated/mean normalized
%points.
X = X - mean(X);
[coeff] = pca(X);
X = X*coeff;
end

