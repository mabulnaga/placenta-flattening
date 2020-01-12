function [ resultImg, Tflat ] = map_intensity_3d( origImg,origMesh, finalMesh, varargin )
%map_intensity_3d: maps the image intensities from the original image space
%to a transformed image space. Pulls intensities using barycentric
%coordinates and trilinear interpolation.
%   Input: origImage: voxel-image with grayscale intensity values
%          origMesh: triangulation mesh
%          finalMesh: triangulation of mesh after transformation
%           varargin: if wanting to force an output image size, place that
%           here.
%   Output: resultImage: voxel-image, same dimension as origImage, that has
%           the grayscale intensity values mapped according to the mesh's
%           transformation. 
%           Tflat: map of where each voxel in the flattened image came from.

D = [1 0 0;0 1 0;0 0 1; -1 -1 -1];
resultImg = zeros(size(origImg)); 

%original image center as your mean!
meshPts = finalMesh.Points ;
range = max(meshPts) - min(meshPts);
%output image size
if(~isempty(varargin))
    imgSize = varargin{1};
else
    imgSize = size(origImg);
    %imgSize = ceil(range + [30,30,10]);
end
%meshPts = meshPts + imgSize/2;
if(sum(min(meshPts)<1)>0) %start at 1, but care about negatives.
    for i = 1:3
        if(min(meshPts(:,i))<1)
            meshPts(:,i) = meshPts(:,i) - min(meshPts(:,i))+5; %prev + 1. Made +5 so its not at the bottom of image
        end
    end
end
%check if mesh goes outside image boundaries.
if(sum((max(meshPts) - imgSize)>-1)>0)
    diff = ceil(max(meshPts) - imgSize +2);
    for i = 1:3
        if(diff(i) > 0)
            imgSize(i) = imgSize(i) + diff(i);
        end
    end
end


finalMesh = triangulation(finalMesh.ConnectivityList,meshPts);
resultImg = zeros(imgSize);
Tflat = zeros([size(resultImg),3]);
%generate coordinates
[xx,yy,zz] = meshgrid(0:size(resultImg,1),0:size(resultImg,2),0:size(resultImg,3));
xxVoxels = xx(1:end-1,1:end-1,1:end-1)+0.5;
yyVoxels = yy(1:end-1,1:end-1,1:end-1)+0.5;
zzVoxels = zz(1:end-1,1:end-1,1:end-1)+0.5;


for i = 1:size(origMesh.ConnectivityList,1)
     xTOrig = origMesh.Points(origMesh.ConnectivityList(i,:),1);
     yTOrig = origMesh.Points(origMesh.ConnectivityList(i,:),2);
     zTOrig = origMesh.Points(origMesh.ConnectivityList(i,:),3);
     
     xT = finalMesh.Points(finalMesh.ConnectivityList(i,:),1);
     yT = finalMesh.Points(finalMesh.ConnectivityList(i,:),2);
     zT = finalMesh.Points(finalMesh.ConnectivityList(i,:),3);
     V= [xT,yT,zT];
     barycentricT = (V'*D);

     maxBox = max([xT,yT,zT],[],1);
     minBox = min([xT,yT,zT],[],1);
     xInside = find(minBox(1)<=xxVoxels & xxVoxels<=maxBox(1));
     yInside = find(minBox(2)<=yyVoxels & yyVoxels<=maxBox(2));
     zInside = find(minBox(3)<=zzVoxels & zzVoxels<=maxBox(3));
     voxels = intersect(intersect(xInside,yInside),zInside);
     voxelCoords = [xxVoxels(voxels)';yyVoxels(voxels)';zzVoxels(voxels)'];

     barycentricCoords = barycentricT\(voxelCoords-repmat([xT(4);yT(4);zT(4)],[1,size(voxelCoords,2)]));
     barycentricCoords = [barycentricCoords;1-sum(barycentricCoords)];
     ind = min(barycentricCoords>=0) & abs(sum(barycentricCoords,1)-1)<=10^-5;
     
     xNew = xT'*barycentricCoords(:,ind);
     yNew = yT'*barycentricCoords(:,ind);
     zNew = zT'*barycentricCoords(:,ind);
     %coord = [round(xNew)',round(yNew)', round(zNew)'];
     %this step is not really necessary. it is the same as
     %voxelCoords(:,ind), then adding 0.5 to get to an integer value for
     %the resulting coorindates.
     coord = [ceil(xNew)', ceil(yNew)', ceil(zNew)'];
     
     xPull = xTOrig'*barycentricCoords(:,ind);
     yPull = yTOrig'*barycentricCoords(:,ind);
     zPull = zTOrig'*barycentricCoords(:,ind);
     
%      for x = xPull
%          xq(voxelCoords??) = x;
%      end
     
     %coordPull = [round(xPull)',round(yPull)',round(zPull)'];
     coordPull = [ceil(xPull)', ceil(yPull)', ceil(zPull)'];

     resultImg(sub2ind(size(resultImg),coord(:,1),coord(:,2),coord(:,3))) = interp3(origImg,yPull,xPull, zPull)';     
     Tflat(sub2ind(size(Tflat),coord(:,1),coord(:,2),coord(:,3),ones(size(coord(:,3))))) = coordPull(:,1);
     Tflat(sub2ind(size(Tflat),coord(:,1),coord(:,2),coord(:,3),2*ones(size(coord(:,3))))) = coordPull(:,2);
     Tflat(sub2ind(size(Tflat),coord(:,1),coord(:,2),coord(:,3),3*ones(size(coord(:,3))))) = coordPull(:,3);

     %old looped code
%      for ii=1:size(coord,1)
%          try
%          %resultImg(coord(ii,1),coord(ii,2),coord(ii,3)) = origImg(coordPull(ii,1), coordPull(ii,2), coordPull(ii,3));
%          Tflat(coord(ii,1),coord(ii,2),coord(ii,3),:) = [coordPull(ii,1),coordPull(ii,2),coordPull(ii,3)];
% 
%          catch
%              1
%          end
%         %         if(origImg(coordPull(ii,1)+1,coordPull(ii,2)+1,coordPull(ii,3)+1))
% %           1  
% %         end
%      end    
end


end

