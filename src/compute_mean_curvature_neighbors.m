function [W] = compute_mean_curvature_neighbors(X, A, T, lambda)
%compute the euclidean distance between a vertex and all of its neighbors.
%Outputs a sparse matrix., exp(-||x_i-x_j||^2/(2\sigma^2))
%X: vertex coordinates
%A: sparse adjacency matrix
%sig: neighborhood weighting
W = zeros(size(A));
[meanCurve,vertexNormals] = compute_mean_curvature(T, X);
    for i = 1:length(A)
        elms = find(A(i,:));
         for ii = 1:length(elms)
             W(i,elms(ii)) = exp(lambda*dot(vertexNormals(i,:),vertexNormals(elms(ii),:))); %metric
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

