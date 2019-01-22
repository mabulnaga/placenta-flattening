function [binMap] = map_boundary_node_index(surfaceTAll,surfaceT)
%Constructs a map between the triangulation from referencing only the nodes
%on the surface (i.e. calling freeBoundary with 2 output arguments), and the
%triangulation that is constructed referencing the nodes in the VOLUME 
%(freeBoundary) with only one output argument.
%Input: surfaceTAll: triangulation with indexes referencing the nodes in
%the original volume
%   surfaceT: triangulation with references to nodes only on the surface
%   (i.e. a new list).
    
binMap = zeros(length(unique(surfaceT)),1);
for i = 1:3
    for j = 1 : length(surfaceT)
        binMap(surfaceT(j,i)) = surfaceTAll(j,i);
    end
end
end

