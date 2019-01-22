function [ vol ] = tet_volume_signed( T )
%Input: tetrahedron, 3x4 matrix. Each row is a vertex, and col 1 is x
%coordinate, col 2 is y coordinate, col 3 is z coordinate. Outputs the
%volume of the tetrahedron.
%Assumes the tetrahedron is 3x4, where each vertex is a column vector.
%vol = abs(dot((T(1,:)-T(4,:)),cross(T(2,:)-T(4,:),...
  %  T(3,:)-T(4,:))))/6;
  D = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
  a = ones(4,1);
  vol = 1/6*(det([T',a]));
 %vol = norm(T*D);
    vol = gather(vol);
end

