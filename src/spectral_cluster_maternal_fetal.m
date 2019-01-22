function [IDX] = spectral_cluster_maternal_fetal(W)
%function to cluster mesh into two, using spectral clustering and cosine
%similarity metric.
%Input: W is a similarity metric (sparse).
%Output: IDX, which is 1 or 2, indicating indices for either cluster.
%Uses the second eigen function to split into two nodal domain.
eigFuncNum = 2;
[row, ~, v] = find(W);
d = accumarray(row, v, []);
D = spdiags(d(:),0,length(d),length(d));
D = spfun(@(x) x.^(-1/2),D);
%compute laplacian.
L=W;
L = D*L*D;
L = speye(length(L))-L;
[V,D] = eigs(L,eigFuncNum,'sa');
V = V(:,eigFuncNum);
%split into two clusters
IDX = zeros(size(V));
IDX(V>0) = 2; IDX(V<=0) = 1;
end

