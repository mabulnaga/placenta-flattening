function [node,elem,face] = create_mesh(img, varargin)
%Function to generate a mesh from a binary segmentation.
%See the iso2mesh documentation for additional details.
%Inputs:
%       img: binary image
%       varargin: meshing parameters. Default: 2, 25.
%Outputs:
%       node: extracted nodes
%       elem: connectivity list
%       face: set of faces
if(nargin == 2)
    [node,elem,face] = v2m(img,0.5,varargin{1},25); 
elseif(nargin == 3)
    [node,elem,face] = v2m(img,0.5,varargin{1},varargin{2});
else
    [node,elem,face] = v2m(img,0.5,2,25); 
    if(length(elem) > 25000)
        %catch to make sure the mesh is not too coarse.
        [node,elem,face] = v2m(img,0.5,3,25); 
    end
end
end