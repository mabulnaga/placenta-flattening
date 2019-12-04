function [ distances ] = compute_geodesic_distance_approximate( X, surfaceT, nodes )
%computes the approximate geodesic distance to nodes adjacent
distances = zeros(length(nodes));
DT = triangulation(surfaceT, X);
E = edges(DT);
edge_len = abs(DT.Points(E(:,1)) - DT.Points(E(:,2)));
h = mean(edge_len);
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
    distances(:,i) = graphshortestpath(weightedA,i)';
end

end
