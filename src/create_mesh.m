function [node,elem,face] = create_mesh(img)
[node,elem,face] = v2m(img,0.5,2,25); %base: 2 (3rd parameter) prev 30
%plotmesh(node,face);
%view(-20,25)
%Modifications: load just image, create mesh, save image to folder. Try
%different parameter settings, build histogram of #tets. How to evaluate?
%maybe convert back to voxel and compare to original??
end