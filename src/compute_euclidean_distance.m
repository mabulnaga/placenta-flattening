function [ distances ] = compute_euclidean_distance( X, surfaceT, nodes )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 distances = zeros(length(nodes));
 DT = triangulation(surfaceT, X);
 E = edges(DT);
 edge_len = abs(DT.Points(E(:,1)) - DT.Points(E(:,2)));
 h = mean(edge_len);
 %X = X(nodes,:);
 for i = 1:length(nodes)
    l = X;
    l(i,:) = [];
    d = vecNorm(l - X(i,:), 2);
    if(i== 1)
        d = [0;d];
    elseif(i <length(nodes))
        d = [d(1:i-1,:); [0]; d(i:end,:)];
    else
        d = [d; [0]];
    end
    distances(:,i) = d;
 end

end

