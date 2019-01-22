function [IDX, nNodes] = find_rim_neighbor_distance(A, P, clusterIDX, neighborDist, varargin)
%Given a cluster of the placenta (2 cluster), finds the rim based on a
%neighborhood distance of multiple nodes. Assumes that if there is a change
%in cluster membership within node i's neighborhood, then node i will be
%assigned to the rim.
%Inputs:
%   A: adjacency matrix (sparse) AKA 1-ring neighborhood
%   P: nvx3: list of vertex coordinate positions
%   clusterIDX: matrix containing 1 and 2, indicating clutser of fetal or
%   maternal side. If it contains 0 or 1 instead, it will convert it.
%   neighborDist: the distance threshold for finding the neighbors of a
%   node. I.e., all nodes with a shortest path < neighborDist are
%   considered in the neighborhood when it comes to finding the rim.
%Outputs:
%    IDX: per-node indices for rim, fetal, and maternal side. Assigns the
%    label '3' to the rim.
%    nNodes: percentage of surface nodes on the rim.

D = compute_euclidean_distance_neighbors(P,A,[],0);
G = graph(D);
d = distances(G);

neighborDistO = neighborDist;
addDist = 0;
%check if input index list has 2 clusters
nEl = length(unique(clusterIDX));
if(nEl ~=2)
    error('more than 2 clusters inputted')
else %check if 0 or 1 inidices.
   i = find(nEl ==0);
   clusterIDX(i) = 2;
end
%Here, compute the rim, and then run some checks.
computeRim = 1;
while(computeRim)
   Adist = zeros(size(A));
   for i = 1:length(A)
       ids = find(d(i,2:end)<=neighborDist);
       ids = ids+1;
       Adist(i,ids) = 1;
   end
    IDX = clusterIDX;
    IDX(find_rim_from_2_cluster(Adist,clusterIDX)==1) = 3; 
    %trisurf(T,P(:,1),P(:,2),P(:,3),IDX);
    %trisurf(T,P(:,1),P(:,2),P(:,3),clusterIDX);
    %now, have the rim. Run some checks. Mainly, that have at least 30% of
    %nodes on the rim (maybe remove later), and that we dont cut the fetal
    %or maternal sides.
    nNodes = sum(IDX==3)/length(IDX);
    binA = conncomp(graph(adjacency_matrix_per_label(A,IDX,1)));
    binB= conncomp(graph(adjacency_matrix_per_label(A,IDX,2)));
    binRim = conncomp(graph(adjacency_matrix_per_label(A,IDX,3)));
    numConA = length(unique(binA(IDX==1)));
    numConB = length(unique(binB(IDX==2)));
    numConRim = length(unique(binRim(IDX==3)));

    if(numConA ~=1 || numConB ~=1)
        computeRim = 1;
        if(addDist == 0)
            neighborDist = neighborDist - 0.5;
        else
            neighborDist = neighborDist + 0.5;
        end
        warning('rim breaks up fetal or maternal side');
%     elseif(numConRim ~=1)
%         computeRim=1;
%         neighborDist = neighborDist + 0.5;
%         warning('rim broken up (non conitnuous');
    elseif(nNodes < 0.30) %probably should remove this condition!
        computeRim = 1;
        neighborDist = neighborDist+0.25;
        warning('too few nodes <threshold %');
    else
        computeRim = 0;
    end
    if(neighborDist <= 0 )
       neighborDist = neighborDistO;
       addDist = 1;
    end
end
end

