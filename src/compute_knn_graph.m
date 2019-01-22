function [A] = compute_knn_graph(A,k)
%expands the adjacency matrix A to its k-nearest neighbors. if k=1, simply
%returns the adjacency matrix.
%inputs: A: adjacency matrix
%k: number of neighbors to expand by.
%To-Do: Make this FASTER and use less loops. Also upload to file exchange.
    B = A;
for i = 1: k-1
    B = B + B*B';
end
A = B;
A(A>1) = 1;
A = A - speye(length(A));
end

