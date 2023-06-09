function [resultImg, Xmap] = pullback_image(T,Xs,Xmap,img_original,interp)
%Pushes forward a map.
%Inputs: T: Ntx4 tetrahedralization
%       Xs: Nx3 source vertices
%       Xmap: Nx3 mapped vertices
%       img_original: 3D image where T, Xs came from
%       We want to create a new image based on Xmap, pulling intensities
%       from Xs back.
%       interp: 'linear', 'nearest', etc...
%Outputs: img size of img_original with pulled back values.
if(nargin == 4)
    interp ='linear';
end
% create a voxel grid same size as img original
imgSize = size(img_original);
% add it to the midpoint...
if(any(min(Xmap)<1))
    Xmap = Xmap + round(imgSize/2);
end
% make sure you don't have any negative values
if(any(min(Xmap)<1))
    round_coord = ceil(abs(min(Xmap)))+1;
    l = min(Xmap) < 1;
    Xmap(:,l) = Xmap(:,l) + round_coord(l);
end
% check if coordinates exceed image bound
rem = imgSize - max(Xmap);
if(any(rem < 1))
    l = rem < 1;
    imgSize(l) = imgSize(l) + ceil(abs(rem(l)))+1;
end
resultImg = zeros(imgSize);
%[xx,yy,zz] = meshgrid(1:size(resultImg,1),1:size(resultImg,2),1:size(resultImg,3));
[xx,yy,zz] = meshgrid(1:imgSize(1),1:imgSize(2),1:imgSize(3));

xxVoxels = xx;
yyVoxels = yy;
zzVoxels = zz;

coord = round([xxVoxels(:), yyVoxels(:), zzVoxels(:)]);

mappedMesh = triangulation(T,Xmap);
[ID,B] = pointLocation(mappedMesh,coord);
%map the voxels forward. Each element of ID says which tet it came from,
%and B is the barycentric coordinates.
voxel_inside = find(~isnan(ID));
xNew = 0;
for i = 1 : 4
    xNew = xNew + Xs(T(ID(voxel_inside),i),:).*B(voxel_inside,i);
end
mappedVoxels = xNew;
resultImg(sub2ind(size(resultImg),coord(voxel_inside,1),coord(voxel_inside,2),coord(voxel_inside,3))) = interp3(img_original,mappedVoxels(:,2),mappedVoxels(:,1),mappedVoxels(:,3),interp,-1)';