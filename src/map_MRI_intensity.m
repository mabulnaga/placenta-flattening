function [mappedImage] = map_MRI_intensity(startVolume, mappedVolume, mriImage)
%Maps MRI signal intensities to the flattened space.
%Inputs:
%       startVolume: Ntx4 triangulation of mesh representing placenta
%           segmentation in source space.
%       mappedVolume: Ntx4 triangulation of mesh representing placenta
%           after flattening
%       mriImage: 3D or 4D matrix containing MRI signals in the source
%       space.
%Returns:
%       mappedImage: matrix containing MRI signals in the flattened space.


dims = size(mriImage);
T = startVolume.ConnectivityList;
Xs = startVolume.Points;
Xmap = mappedVolume.Points;
if length(dims)==3 || dims(4)==1 %image is 3D or 4th dimension is singular
    [mappedImage, ~] = pullback_image(T,Xs,Xmap,double(mriImage),'linear');
elseif dims(4)>1 %image is 4D - map the intesity values by looping over volumes
    mappedImage = [];
    for i=1:dims(4) 
        [mappedImageTmp, ~] = pullback_image(T,Xs,Xmap,double(mriImage(:,:,:,i)),'linear');
        mappedImage = cat(4,mappedImage,mappedImageTmp);
    end         
end

end

