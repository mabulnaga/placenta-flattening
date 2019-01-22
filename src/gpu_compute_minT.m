
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
  
  %Going to try log-sum exp. Remove stuff with Z later.
%   d1 = -(a13.*a22.*a31); d2 = a12.*a23.*a31; d3 = a13.*a21.*a32; d4 = -a11.*a23.*a32; d5 = -a12.*a21.*a33; d6 = a11.*a22.*a33;
%   d1 = reshape(d1,[size(d1,3),1]);
%   d2 = reshape(d2,[size(d2,3),1]);
%   d3 = reshape(d3,[size(d3,3),1]);
%   d4 = reshape(d4,[size(d4,3),1]);
%   d5 = reshape(d5,[size(d5,3),1]);
%   d6 = reshape(d6,[size(d6,3),1]);
%   v = [d1,d2,d3,d4,d5,d6];
%   max_v = max(v,[],2);
%   min_v = min(v,[],2);
%   maxes = [max_v,abs(min_v)];
%   vals = [max_v, min_v];
%   [x,ind] = max(maxes,[],2);
%   zz = length(ind).*(ind-1);
%   l = 1:length(ind);
%   l = l';
%   zz = zz+l;
%   constant(:,1) = vals(zz);
%   e = log(sum(exp(v-constant),2))+constant;
  
  d = -(a13.*a22.*a31) + a12.*a23.*a31 + a13.*a21.*a32 - a11.*a23.*a32 - a12.*a21.*a33 + a11.*a22.*a33;
  c = (-(a23.*a32.*b11) + a22.*a33.*b11 + a23.*a31.*b12 - a21.*a33.*b12 - a22.*a31.*b13 + a21.*a32.*b13 + a13.*a32.*b21 - a12.*a33.*b21 - a13.*a31.*b22 + a11.*a33.*b22 + a12.*a31.*b23 - a11.*a32.*b23 - a13.*a22.*b31 + a12.*a23.*b31 + a13.*a21.*b32 - a11.*a23.*b32 - a12.*a21.*b33 + a11.*a22.*b33);
  a = (-(b13.*b22.*b31) + b12.*b23.*b31 + b13.*b21.*b32 - b11.*b23.*b32 - b12.*b21.*b33 + b11.*b22.*b33);
  b = (-(a33.*b12.*b21) + a32.*b13.*b21 + a33.*b11.*b22 - a31.*b13.*b22 - a32.*b11.*b23 + a31.*b12.*b23 + a23.*b12.*b31 - a22.*b13.*b31 - a13.*b22.*b31 + a12.*b23.*b31 - a23.*b11.*b32 + a21.*b13.*b32 + a13.*b21.*b32 - a11.*b23.*b32 + a22.*b11.*b33 - a21.*b12.*b33 - a12.*b21.*b33 + a11.*b22.*b33);
  %vol = abs(det(T*D));
  %vol = norm(T*D);
  %vol = gather(vol);

  a = -1 * a;
  c = -1 * c;
  %had these because I derived the wrong thing in mathematica. 
%multiply all coefficients by the global minimum. Didn't work very well.
%   z = [min(abs(a(abs(a)>0))), min(abs(b(abs(b)>0))), min(abs(c(abs(c)>0))), min(abs(d(abs(d)>0)))];
%   z = min(z);

%multiply each coefficient by the minimum, to normalize it. Had errors.
%Not sure if I need this, but will add back later if I do.
%   z = [a b c d];
%   z(~z) = Inf;
%   z = min(abs(z));
%   z = pagefun(@inv,z);
% % %   z=1e10;
%   a = pagefun(@mtimes,a,z);
%   b = pagefun(@mtimes,b,z);
%   c = pagefun(@mtimes,c,z);
%   d = pagefun(@mtimes,d,z);

%simple solution, one fixed coefficient. No real erasoning, but worked well
%emprically.
%    z = 1e-10;
%    a=a./z;
%    b=b./z;
%    c=c./z;
%    d=d./z;

  disc = 18*a.*b.*c.*d - 4*b.^3.*d +b.^2.*c.^2 -4*a.*c.^3 - 27*a.^2.*d.^2;
  delta0 = b.^2 - 3*a.*c;
  delta1 = 2*b.^3 -9*a.*b.*c+27*a.^2.*d;
%   multiroots = zeros(length(disc),1);
%   distroots = zeros(length(disc),1);
%   complexroots = zeros(length(disc),1);
%   singleroots = zeros(length(disc),1);
%   
%   multiroots(find(disc==0)) = 1;
%   distroots(find(disc>1e-10)) = 1;
%   complexroots(find(disc<0)) = 1;
  
  roots1 = Inf*ones(length(disc),3);
  roots1 = gpuArray(roots1);
  roots2 = roots1;
  roots3 = roots1;
  delta0 = gather(delta0);
  delta1 = gather(delta1);
  
%   singleroots(find(delta0<1e-80 & delta0 >0))=1;
%   tworoots = not(singleroots);
%   singleroots = multiroots & singleroots;
%   tworoots = multiroots & tworoots;
  
  
%   roots1(singleroots) = gather(-b(singleroots)/3.*(a(singleroots).^(-1)));
%   roots2(singleroots) = -Inf;
%   roots3(singleroots) = -Inf;
%   
%   roots1(tworoots) = gather(1/2*(9*a(tworoots).*d(tworoots)-b(tworoots).*c(tworoots)).*delta0(tworoots).^(-1));
%   roots2(tworoots) = gather((4*a(tworoots).*b(tworoots).*c(tworoots)-9*a(tworoots).^2.*d(tworoots)-b(tworoots).^3).*a(tworoots).^(-1).*delta0(tworoots).^(-1));
%   roots3(tworoots) = -Inf;
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
  r= gather(min(rootsall));
  rmin = r;
  digits(digitsOld);
  %% 
%    x = find(rootsall==min(rootsall));
%    rminC = rootsall(x);
%    [row, col] = ind2sub(size(rootsall),x); 
%    x =row;
%    aa = a(x);
%    bb = b(x);
%    cc = c(x);
%    dd = d(x);
%    C1 = C(x);
%    delta01 = delta0(x);
%    delta11 = delta1(x);
%    disc1 = disc(x);
% %   
% % %Old method with matlab's roots
% roots11=Inf*ones(length(disc),3);
% roots11= gpuArray(roots11);
% roots111 = roots11;
% for ii=1:length(a)
%     r = roots([a(ii),b(ii),c(ii),d(ii)]);
%     try
%         roots111(ii,1:length(r)) = r;
%     catch
%         l =1;
%     end
%     %r = r(imag(r)==0);
%     %r = r(r>0);
%     r(imag(r)~=0) = Inf;
%     r(r<0) = Inf;
%     if(length(r)>=1)
%         roots11(ii,1:length(r))=r;
%     end
% end
% r = gather(min(min(roots11)));
% x = find(roots11==min(min(roots11)));
% [row, col] = ind2sub(size(roots11),x); 
% rC = roots111(row,col);
%   
% if(isinf(r))
%     r=10;
% end
% if(abs(r-rmin)>1e-5)
% fprintf ('Roots min: %d, quick min: %d, a:%1.2d, b:%1.2d, c:%1.2d, d:%1.2d, C:%1.2d, disc: %1.2d, delta0: %1.2d, delta1: %1.2d\n', r, rmin,aa,bb,cc,dd, C1, disc1, delta01, delta11)
% v = gather([rmin-r, rC, r, rminC, rmin, aa, bb, cc,dd, C1, disc1, delta01, delta11]);
% dlmwrite('Myfile.csv',v,'-append');
% end

%Checks
% diffs = zeros(2341,3);
% for jj=1:2341
%     diffs(jj,1) =gather( roots1(jj,1)-roots11(jj));
%     diffs(jj,2) = gather(roots1(jj,2)-roots2(jj));
%     diffs(jj,3) = gather(roots1(jj,3) - roots3(jj));
% end
%Cardano's Method
% c1=c;c0=d;c2=b;c1=a; 
% p = (3*c1-c2.^2)/3;
% q = (9*c1.*c2-27*c0-2*c2.^3)/27;
% Q=(1/3)*p;
% R=(1/2)*q;
% D=Q.^3+R.^2;
% k=find(D>=0);
% S=zeros(length(D),1);
% T=zeros(length(D),1);
% S(k)=nthroot(gather(R(k)+sqrt(D(k))),3);
% T(k)=nthroot(gather(R(k)-sqrt(D(k))),3);
% root=zeros(length(D),1);
% root(k,1) = -(1/3)*c2+(S(k)+T(k));
% root(k,2) = -(1/3)*c2-(1/2)*(S(k)+T(k))+(1/2)*i*sqrt(3)*(S(k)-T(k));
% root(k,3) = -(1/3)*c2-(1/2)*(S(k)+T(k))-(1/2)*i*sqrt(3)*(S(k)-T(k));
% k=find(D<0);
% phi=zeros(length(D),3);
% phi(k)=acos(R(k)./sqrt(-Q(k).^3));
% root(k,1) = 2*sqrt(-Q(k)).*cos(phi(k)/3)-(1/3)*c2;
% root(k,2) = 2*sqrt(-Q(k)).*cos((phi(k)+2*pi)/3)-(1/3)*c2;
% root(k,3) = 2*sqrt(-Q(k)).*cos((phi(k)+4*pi)/3)-(1/3)*c2;
  
  
end
