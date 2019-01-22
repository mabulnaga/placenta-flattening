function [A] = adjacency_matrix_per_label(A, labels, labelID)
%ADJACENCY_MATRIX_PER_LABEL outputs a square adjacency matrix, including
%only the nodes with a certain label ID.
%Inputs:
%   A: sparse adjacency matrix
%   labels: vector size nv x 1, containing the set of labels corresponding
%   to each node
%   labelID: scalar, includes the one label ID we want to keep
for i = 1: length(A)
    if(labels(i)~=labelID)
        [k,j] = find(A(i,:));
        A(i,:)= 0;
        A(j,i) = 0;
    end
    
end

