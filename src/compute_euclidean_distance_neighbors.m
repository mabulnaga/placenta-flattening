function [W] = compute_euclidean_distance_neighbors(X, A, sig, varargin)
%compute the euclidean distance between a vertex and all of its neighbors.
%Outputs a sparse matrix., exp(-||x_i-x_j||^2/(2\sigma^2))
%X: vertex coordinates
%A: sparse adjacency matrix
%sig: neighborhood weighting
useExp = 0;
W = zeros(size(A));
    for i = 1:length(A)
        elms = find(A(i,:));
        for ii = 1:length(elms)
            W(i,elms(ii)) = vecNorm(X(i,:) - X(elms(ii),:),2); %metric
        end
    end
if(nargin == 4)
    useExp = varargin{1};
end
W = sparse(W);
if(useExp)
    W = spfun(@(x) -1.*x.^2./(2*sig^2),W);
    W = spfun(@exp,W);
end
end

