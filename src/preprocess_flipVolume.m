function [ T ] = preprocess_flipVolume( T, X )
%PREPROCESS_FLIPVOLUME Summary of this function goes here
%   Detailed explanation goes here
%Input a tet, output a tet with all positive signed volumes.
try
    canUseGPU = parallel.gpu.GPUDevice.isAvailable;
catch 
    canUseGPU = false;
end
if(canUseGPU)
    X0 = gpu_generate_tets(T, X);
else
    X0 = cpu_generate_tets(T, X);
end
Xvol = arrayfun(@(x) tet_volume_signed(X0(:,:,x)), 1:size(X0,3));
indices = find(Xvol<0);
for i = 1: length(indices)
    l = T(indices(i),2);
    l2 = T(indices(i),4);
    T(indices(i),2) = l2;
    T(indices(i),4) = l;
end

end

