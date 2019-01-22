function [ distortion ] =compute_distortion( X,Xbar,D)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    %J = X'*D\(Xbar'*D);
    %J = (D*X')\((D*Xbar'));
    
   % D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];
    %J = X'*D*(Xbar'*D)^(-1);
%      X = [X'; 1 1 1 1];
%     Xbar = [Xbar'; 1 1 1 1];
%     c = [1 0 0 0;0 1 0 0;0 0 1 0];
%     b= [1 0 0;0 1 0;0 0 1;0 0 0];
%     J = c*X*inv(Xbar)*b;
    distortion = norm(J,'fro')^2 + norm(inv(J),'fro')^2;
end

