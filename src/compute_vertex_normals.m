function [ vertexNormals] = compute_vertex_normals( T, boundaryVertices, faceNormals )
%input: T: list of triangles Nx3, where the indices are the same as in
%boundary Vertices
%   boundaryVertices: unique list of nodes, size Npx1
%   areas: surface areas, size 1xN
%   Detailed explanation goes here
%output: voronoi area per node, size 1xNp
L = zeros(length(T),length(boundaryVertices));
L = logical(L);
vertexNormals = zeros(length(boundaryVertices),3);
for i = 1: length(boundaryVertices)
    v1 = ismember(T(:,1),boundaryVertices(i));
    v2 = ismember(T(:,2),boundaryVertices(i));
    v3 = ismember(T(:,3),boundaryVertices(i));
    L(:,i) = v1 | v2 | v3;
    %voronoi_areas(i) = accumarray(find(L(:,i))',areas');
    %vertexNormals(i) = sum(areas(L(:,i)).*faceNormals(L(:,i),:));
    %vertexNormals(i) = sum(repmat(areas(L(:,i))',1,3).*faceNormals(L(:,i),:),1)/sum(areas(L(:,i)));
    vertexNormals(i,:) = sum(faceNormals(L(:,i),:),1);
end
%vertexNormals = vertexNormals / 3;
