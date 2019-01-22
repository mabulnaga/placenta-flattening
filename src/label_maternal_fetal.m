function [maternalNodes, fetalNodes] = label_maternal_fetal(T, P, northNodes, southNodes, varargin)
%Given the north and south side segmentations (rim or no rim), this
%function assigns the label maternal to the south nodes, thereby allowing
%the algorithm to always have the fetal side on top.
%Inputs:
%       T: surface triangulation
%       P: surface points
%       northNodes: indices of points with the label 'North'
%       southNodes: indices of points with the label 'South'
%       varargin: used for debugging. For now, varargin{1} is the
%       voronoiAreas, for a weighted mean (depricated).
%Outputs:
%       maternalNodes: indices of nodes with the label maternal
%       fetalNodes: indices of nodes with the label fetal


K = convhull(P(:,1),P(:,2),P(:,3),'simplify',true);
nNorth = sum(ismember(northNodes,K));
nSouth = sum(ismember(southNodes,K));
%the maternal side will have more nodes on the hull.
if (nNorth > nSouth)
    maternalNodes = northNodes;
    fetalNodes = southNodes;
else
    maternalNodes = southNodes;
    fetalNodes = northNodes;
end

%% Depricated methods
% Clustering based on average mean curvature. This was too mesh dependent.
% meanCurvature = compute_mean_curvature(T,P);
% if(~isempty(varargin))
%     meanCurvature = meanCurvature .*varargin{1};
%     m1 = sum(meanCurvature(northNodes));
%     m2 = sum(meanCurvature(southNodes));
% else
%     m1 = mean(meanCurvature(northNodes));
%     m2 = mean(meanCurvature(southNodes));
% end
% 
% 
% 
% if (m1 > m2)
%     maternalNodes = northNodes;
%     fetalNodes = southNodes;
% else
%     maternalNodes = southNodes;
%     fetalNodes = northNodes;
% end

%Old: testing using the convex hull:
% K = convhull(P(:,1),P(:,2),P(:,3));
% KN = convhull(P(northPts,1),P(northPts,2),P(northPts,3));
% KS = convhull(P(southPts,1),P(southPts,2),P(southPts,3));
% trisurf(T,P(:,1),P(:,2),P(:,3))
% hold on
% plot3(P(K(:,1),1),P(K(:,2),2),P(K(:,3),3),'ro')
% hullPts = [P(K(:,1),1),P(K(:,2),2),P(K(:,3))];
% hullPtsN = P(K(:,1),1),P(K(:,2),2),P(K(:,3))
% northDists = zeros(1,length(hullPts));
% southDists = zeros(1,length(hullPts));
% for i = 1:length(northNodes)
%     p = repmat(P(northNodes(i),:),length(hullPts),1);
%     d = vecNorm(p-hullPts,2,2);
%     northDists(i) = min(d);
% end
% for i = 1:length(southNodes)
%     p = repmat(P(southNodes(i),:),length(hullPts),1);
%     d = vecNorm(p-hullPts,2,2);
%     southDists(i) = min(d);
% end
% %Want: fetal side to have longer distance to hull.
% m1Conv = mean(northDists);
% m2Conv = mean(southDists);


end

