function [mappedImage] = main(grayImage, segmentationImage, varargin)
%Main function to call placenta parameterization code. Takes as input a
%gray MRI image (ex: BOLD, HASTE, etc...), either as a 3D image file (a
%matrix), or as a NIFTI file, and the corresponding placenta segmentation
%and maps it to the flattened space. The mapping produces a flattened mesh
%and a transformed grayscale image. This script will also produce a blank NIFTI
%image with the transformed output (binary and grayscale).
%Note that the output will only be a NIFTI file if the inputs are NIFTI
%files.
%Inputs:
%       grayImage: MRI image (HASTE, BOLD, etc), either as a 3D matrix, or
%       as a NIFTI file. If Nifti file, input as a string to the full path.
%       segmentationImage: corresponding binary image where the placenta is
%       segmented. Either nifti or image file (3D).
%       varargin:
%               1. savePath. Specify a full path to save the output to. If
%               left empty or unspecified, by default will create a folder outputs/ and save there.
%Outputs:
%       mappedImage: the mapped MRI image to the flattened space (matrix).
%       Output is a NIFTI file if the inputs were NIFTI files.

%Default algorithm parameters
lambda = 1;
rho = 1/2;
saveNifti = 0;
%Load in images and check if they are NIFTI files


if(ischar(grayImage))
    try
        grayImageNifti = loadNii(grayImage);
        grayImage = grayImageNifti.img;
        saveNifti = 1;
    catch
        error('invalid nifti file specified');
    end
end
if(ischar(segmentationImage))
    try
        segmentationImageNifti = loadNii(segmentationImage);
        segmentationImage = segmentationImageNifti.img;
        segmentationImage = preprocess_seg_holes(segmentationImage);
    catch
        error('invalid nifti file specified');
    end
else
    segmentationImage = preprocess_seg_holes(segmentationImage);
end

%Check if the image mask is binary
if(length(unique(segmentationImage))>2)
    error('must use a binary image mask')
end

%create a directory to save output files to
currentDir = mfilename('fullpath');
currentDir = fileparts(currentDir); currentDir = fileparts(currentDir);
if(~isempty(varargin))
    savePath = varargin{1};
    try
        mkdir(savePath)
    catch
        savePath = [currentDir,'/output/'];
        if(~exist(savePath,'dir'))
           mkdir(savePath); 
        end
        savePath = [savePath,datestr(now,'dd-mm-yyyy-HH-MM-SS'),'/'];
        if(~exist(savePath,'dir'))
            mkdir(savePath);
        end
    end
else
    savePath = [currentDir,'/output/'];
    if(~exist(savePath,'dir'))
        mkdir(savePath);
    end
    savePath = [savePath,datestr(now,'dd-mm-yyyy-HH-MM-SS'),'/'];
    if(~exist(savePath,'dir'))
        mkdir(savePath);
    end
end

%mapping to the flattened space
[startVolume, ~, mappedVolume] = flatten_algorithm(segmentationImage, lambda, rho);

%saving the mapped meshes
fprintf('Flattening complete, saving meshes to %s \n', savePath);
save([savePath,'/startVolume'],'startVolume','-v7.3');
save([savePath,'/mappedVolume'],'mappedVolume','-v7.3');
fprintf('Mapping MRI image to flattened space...\n');

%map the MRI intensity values to the flattened space
[mappedImage, Tflat] = map_intensity_3d(grayImage, startVolume, mappedVolume);

%save the mapped images as NIFTI files if NIFTI files were input to the
%algorithm. Always saves outputs as matlab matrices.
if(saveNifti==1)
    fprintf('saving mapped image as a NIFTI file \n');
    grayImageNifti = update_nifti(grayImageNifti, mappedImage);
    segMapped = logical(mappedImage);
    segImageNifti = update_nifti(grayImageNifti, segMapped);
    saveNii(grayImageNifti, [savePath,'flat-MRI.nii']);
    saveNii(segImageNifti,[savePath,'flat-segmentation.nii']);
    save([savePath,'flat-MRI.mat'],'mappedImage','-v7.3');
    save([savePath,'flat-segmentation.mat'],'segMapped','-v7.3');
else
    fprintf('saving mapped image as a matlab matrix file');
    save([savePath,'flat-MRI.mat'],'mappedImage','-v7.3');
    segMapped = logical(mappedImage);
    save([savePath,'flat-segmentation.mat'],'segMapped','-v7.3');
    niftiwrite(mappedImage,[savePath,'flat-MRI.nii']);
    niftiwrite(logical(mappedImage),[savePath,'flat-segmentation.nii']);
end

end

