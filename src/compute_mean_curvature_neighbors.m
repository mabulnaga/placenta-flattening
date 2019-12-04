function [W] = compute_mean_curvature_neighbors(X, A, T, lambda, varargin)
%compute the euclidean distance between a vertex and all of its neighbors.
%Outputs a sparse matrix., exp(-||x_i-x_j||^2/(2\sigma^2))
%X: vertex coordinates
%A: sparse adjacency matrix
%lambda: neighborhood weighting
%varargin: 1: D: matrix of geodesic distiances
D = ones(size(A));
if(nargin > 4)
   D = varargin{1}; 
end
W = zeros(size(A));
[meanCurve,vertexNormals] = compute_mean_curvature(T, X);
    for i = 1:length(A)
        elms = find(A(i,:));
         for ii = 1:length(elms)
             W(i,elms(ii)) = exp(D(i,elms(ii))*lambda*dot(vertexNormals(i,:),vertexNormals(elms(ii),:))); %metric
         end
%          for ii = 1:length(A)
%              W(i,ii) = dot(vertexNormals(i,:),vertexNormals(ii,:)); %metric
%          end
    end
W = (W + W.')/2;
W = sparse(W);
%W = spfun(@(x) -1.*x.^2./(2*sig^2),W);
%W = spfun(@exp,W);
end

