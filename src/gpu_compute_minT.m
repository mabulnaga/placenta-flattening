function [ r ] = gpu_compute_minT( X,gradX )
%Input: tetrahedron, 4x3 matrix. Each row is a vertex, and col 1 is x
%coordinate, col 2 is y coordinate, col 3 is z coordinate. Outputs the
%volume of the tetrahedron.
%Assumes the tetrahedron is 3x4, where each vertex is a column vector.
%vol = abs(dot((T(1,:)-T(4,:)),cross(T(2,:)-T(4,:),...
  %  T(3,:)-T(4,:))))/6;
  digitsOld = digits(100);
  D = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
  %zeroThresh = eps;
  zeroThresh = 1e-4;
  %LOWERED THRESHOLD TO FIND BUG
  A = pagefun(@mtimes,X,D);
  B = pagefun(@mtimes,gradX,D);
  A = pagefun(@transpose,A);
  B = pagefun(@transpose, B);
  a11 = A(1,1,:);a12 =A(1,2,:);a13 = A(1,3,:);a21 = A(2,1,:);a22 = A(2,2,:); 
  a23= A(2,3,:); a31 = A(3,1,:); a32 = A(3,2,:); a33= A(3,3,:);
  b11 = B(1,1,:);b12 =B(1,2,:);b13 = B(1,3,:);b21 = B(2,1,:);b22 = B(2,2,:); 
  b23= B(2,3,:); b31 = B(3,1,:); b32 = B(3,2,:); b33= B(3,3,:);
  
  
  d = -(a13.*a22.*a31) + a12.*a23.*a31 + a13.*a21.*a32 - a11.*a23.*a32 - a12.*a21.*a33 + a11.*a22.*a33;
  c = (-(a23.*a32.*b11) + a22.*a33.*b11 + a23.*a31.*b12 - a21.*a33.*b12 - a22.*a31.*b13 + a21.*a32.*b13 + a13.*a32.*b21 - a12.*a33.*b21 - a13.*a31.*b22 + a11.*a33.*b22 + a12.*a31.*b23 - a11.*a32.*b23 - a13.*a22.*b31 + a12.*a23.*b31 + a13.*a21.*b32 - a11.*a23.*b32 - a12.*a21.*b33 + a11.*a22.*b33);
  a = (-(b13.*b22.*b31) + b12.*b23.*b31 + b13.*b21.*b32 - b11.*b23.*b32 - b12.*b21.*b33 + b11.*b22.*b33);
  b = (-(a33.*b12.*b21) + a32.*b13.*b21 + a33.*b11.*b22 - a31.*b13.*b22 - a32.*b11.*b23 + a31.*b12.*b23 + a23.*b12.*b31 - a22.*b13.*b31 - a13.*b22.*b31 + a12.*b23.*b31 - a23.*b11.*b32 + a21.*b13.*b32 + a13.*b21.*b32 - a11.*b23.*b32 + a22.*b11.*b33 - a21.*b12.*b33 - a12.*b21.*b33 + a11.*b22.*b33);


  disc = 18*a.*b.*c.*d - 4*b.^3.*d +b.^2.*c.^2 -4*a.*c.^3 - 27*a.^2.*d.^2;
  delta0 = b.^2 - 3*a.*c;
  delta1 = 2*b.^3 -9*a.*b.*c+27*a.^2.*d;

  
  roots1 = Inf*ones(length(disc),3);
  roots1 = gpuArray(roots1);
  roots2 = roots1;
  roots3 = roots1;

  
%   
  C = ((delta1+(delta1.^2-4*delta0.^3).^(1/2))/2).^(1/3);
  zee = -1/2 + 1/2*sqrt(3)*1j;
  roots1 = -1/3*(a.^(-1)).*(b+zee^0*C + 1/(zee^0)*C.^(-1).*delta0);
  roots2 = -1/3*(a.^(-1)).*(b+zee^1*C + 1/(zee^1)*C.^(-1).*delta0);
  roots3 = -1/3*(a.^(-1)).*(b+zee^2*C + 1/(zee^2)*C.^(-1).*delta0);
  rootsall = [roots1(:), roots2(:), roots3(:)];
  %New stuff, extracting real and positive roots only.
  roots1 = roots1(abs(imag(roots1))<zeroThresh);
  roots1 = roots1(roots1>0);
  roots2 = roots2(abs(imag(roots2))<zeroThresh);
  roots2 = roots2(roots2>0);
  roots3 = roots3(abs(imag(roots3))<zeroThresh);
  roots3 = roots3(roots3>0);
  rootsall = [roots1(:); roots2(:); roots3(:)];
  %r= gather(min(rootsall));
  r = min(rootsall);
  rmin = r;
  digits(digitsOld);
  
  
end
