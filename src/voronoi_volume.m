function [ voronoi_volumes ] = voronoi_volume( T, X )
%input: T: list of tets Nx4, where the indices are the same as in
%boundary Vertices
%   X: vertices, size Nvx3. Note that Nv must = max(T).
%   Detailed explanation goes here
%output: voronoi volume per vertex, size 1xNp (NOT NORMALIZED)
vertices = unique(T);
L = zeros(length(T),length(vertices));
L = logical(L);
X0 = cpu_generate_tets(T,X);
Xvol = arrayfun(@(x) tet_volume(X0(:,:,x)), 1:size(X0,3));
voronoi_volumes = zeros(length(X),1);
for i = 1: length(vertices)
    v1 = ismember(T(:,1),vertices(i));
    v2 = ismember(T(:,2),vertices(i));
    v3 = ismember(T(:,3),vertices(i));
    v4 = ismember(T(:,4), vertices(i));
    L(:,i) = v1 | v2 | v3 | v4;
    %voronoi_areas(i) = accumarray(find(L(:,i))',areas');
    voronoi_volumes(vertices(i)) = sum(Xvol(L(:,i)));
end
voronoi_volumes = voronoi_volumes /4;
