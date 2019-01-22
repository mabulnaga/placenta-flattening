function [ area ] = triangle_area( T, X )
%Input: T: set of triangles, with coordinates referencing X
%X: cartesian coordinates of each node
%output: area per triangle
  %  T(3,:)-T(4,:))))/6;
  area = zeros(1,size(T,1));
  for i = 1:size(T,1)
      x1 = X(T(i,1),:);
      x2 = X(T(i,2),:);
      x3 = X(T(i,3),:);
      area(i) = 1/2*norm(cross(x3-x1,x2-x1));
  end

end

