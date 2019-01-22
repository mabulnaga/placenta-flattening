function [outputNifti] = update_nifti(inputNifti, updateImg)
%Function overwrites a NIFTI struct with a new image, taking care to handle
%the sizing and resets the orientation to be at the origin.
%Inputs:
%       inputNifti: input NIFTI struct
%       updateImg: image to ammend the nifti file with
%Outputs:
%       outputNifti: updated NIFTI struct with the new image

outputNifti = inputNifti;
outputNifti.img = updateImg;
qaternBold = outputNifti.hdr.hist.quatern_b;
quaternCold = outputNifti.hdr.hist.quatern_c;
quaternDold = outputNifti.hdr.hist.quatern_d;
outputNifti.hdr.hist.quatern_b = 0;
outputNifti.hdr.hist.quatern_c = 0;
outputNifti.hdr.hist.quatern_d = 0;
outputNifti.hdr.dime.dim(2) = size(updateImg,1);
outputNifti.hdr.dime.dim(3) = size(updateImg,2);
outputNifti.hdr.dime.dim(4) = size(updateImg,3);
end

