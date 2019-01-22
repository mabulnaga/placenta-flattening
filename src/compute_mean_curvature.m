function [meanCurvature,vertexNormals, L, A] = compute_mean_curvature(T,P)
surfaceAreas = triangle_area(T,P);
voronoiAreas = voronoi_area(T,unique(T),P);
L = cotmatrix(P,T);
faceNormals = normals(P,T);
vertexNormals = compute_vertex_normals(T, unique(T), faceNormals);
vertexNormals = vertexNormals./vecNorm(vertexNormals,2,2);
g = 1/2*L*P;
meanCurvature = zeros(length(P),1);
A = spdiags(voronoiAreas, 0, length(P), length(P));
for i = 1:length(P)
    meanCurvature(i) = dot(g(i,:),vertexNormals(i,:))/voronoiAreas(i);
end
end