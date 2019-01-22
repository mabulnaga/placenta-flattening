function [ voronoi_areas ] = voronoi_area( T, boundaryVertices, X )
%input: T: list of triangles Nx3, where the indices are the same as in
%boundary Vertices
%   boundaryVertices: unique list of nodes, size Npx1
%   X: vertices, size Nvx3. Note that Nv must = max(T).
%   Detailed explanation goes here
%output: voronoi area per node, size 1xNp
L = zeros(length(T),length(boundaryVertices));
L = logical(L);
areas = triangle_area(T, X);
voronoi_areas = zeros(length(X),1);
for i = 1: length(boundaryVertices)
    v1 = ismember(T(:,1),boundaryVertices(i));
    v2 = ismember(T(:,2),boundaryVertices(i));
    v3 = ismember(T(:,3),boundaryVertices(i));
    L(:,i) = v1 | v2 | v3;
    %voronoi_areas(i) = accumarray(find(L(:,i))',areas');
    voronoi_areas(boundaryVertices(i)) = sum(areas(L(:,i)));
end
voronoi_areas = voronoi_areas / 3;

