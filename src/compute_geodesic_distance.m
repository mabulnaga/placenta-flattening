function [ distances ] = compute_geodesic_distance( X, surfaceT, nodes )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 distances = zeros(length(nodes));
 DT = triangulation(surfaceT, X);
 E = edges(DT);
 edge_len = abs(DT.Points(E(:,1)) - DT.Points(E(:,2)));
 h = mean(edge_len);
 %X = X(nodes,:);
 A = adjacency_matrix(surfaceT);
 weightedA = A;
 for i = 1:length(A)
     l = find(A(i,:));
     pt = X(i,:);
     pt = repmat(pt,[length(l),1]);
     pt2 = X(l,:);
     norms = vecNorm(pt - pt2,2,2);
    weightedA(i,l) = norms;
 end
 for i = 1:length(nodes)
     if(size(surfaceT,2)>3)
         distances(:,i) = graphshortestpath(weightedA,i)';
     else
        [distances(:,i), S, Q] = perform_fast_marching_mesh(X, surfaceT,nodes(i));
     end
     %distances(:,i) = heat_geodesic(X, surfaceT, nodes(i),h^3);
 end

end

