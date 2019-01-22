function [IDX] = find_rim_from_2_cluster(A, IDX)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lbls = unique(IDX);
IDX2 = repmat(IDX,1,size(A,2));
labels = A.*IDX2';
borderN = zeros(length(IDX),1);
for i = 1: length(labels)
    if(sum(labels(i,:)) ~= length(find(labels(i,:)))*lbls(1) && sum(labels(i,:)) ~=length(find(labels(i,:)))*lbls(2))
        borderN(i) = 1;
    end
end
IDX = borderN;
c=0;
% for j = 1:length(labels)
% if(length(unique(labels(j,:)))>2)
% c = c+1;
% end
% end
end

