function [ tets ] = gpu_generate_tets( T,X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tets = zeros(3,4,size(T,1));
% tets = arrayfun(@(x) create_tet(T(x,:),X), 1:size(T,1),'UniformOutput', false);
% tets = reshape(cell2mat(tets),[3,4,size(T,1)]);
% tets = gpuArray(tets);
% tets2 = tets;
p=1:length(T);
p=p';
tets = gpuArray(zeros(size(T,1)*3,size(T,2)));
tets(:,1) = reshape((X(T(p,1),:))',[],1); 
tets(:,2) = reshape((X(T(p,2),:))',[],1); 
tets(:,3) = reshape((X(T(p,3),:))',[],1); 
tets(:,4) = reshape((X(T(p,4),:))',[],1); 
tets = reshape(tets',4,3,[]);
tets = pagefun(@transpose, tets);
end

