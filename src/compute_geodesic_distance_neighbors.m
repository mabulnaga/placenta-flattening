function [W] = compute_geodesic_distance_neighbors(X, A, Ts, sig)
%compute the euclidean distance between a vertex and all of its neighbors.
%Outputs a sparse matrix., exp(-||x_i-x_j||^2/(2\sigma^2))
%X: vertex coordinates
%A: sparse adjacency matrix
%sig: neighborhood weighting

%Using graph shortest path to compute geodisc distance...
weightedA = A;
%  for i = 1:length(A)
%      l = find(A(i,:));
%      pt = X(i,:);
%      pt = repmat(pt,[length(l),1]);
%      pt2 = X(l,:);
%      norms = vecNorm(pt - pt2,2,2);
%     weightedA(i,l) = norms;
%  end
 W = zeros(size(A));
    for i = 1:length(A)
        elms = find(A(i,:));
        d = perform_fast_marching_mesh(X,Ts,i);
        W(i,elms) = d(elms);
        %for ii = 1:length(elms)
            %W(:,i) = graphshortestpath(weightedA,i)'; %metric
        %end
    end
W = sparse(W);
%W = spfun(@(x) -1.*x.^2./(2*sig^2),W);
%W = spfun(@exp,W);
end

