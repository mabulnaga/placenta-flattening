function grad = distortion_grad(X,Xbar)
    D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
    Xbar = Xbar';
    O = Xbar*D;
    Ebar = 0;
    X = X';
    %grad = X*D*(inv(Ebar).'*inv(Ebar)+inv(Ebar)*inv(Ebar.'))...
   %     -X*D*inv(D.'*X.'*X*D)*(2*Ebar.'*Ebar)*inv(D.'*X.'*X*D);
%     grad = 2*X*D*(inv(Ebar)*inv(Ebar).')...
%         -X*D*inv(D.'*X.'*X*D)*(2*Ebar.'*Ebar)*inv(D.'*X.'*X*D);
    %grad = X*(D.'*D*inv(O)*inv(O).')+X*(D.'*D*inv(O)*inv(O).')-...
    %    X*inv(X.'*X)*(inv(D)*inv(D).'*O.'*O+inv(inv(D)*inv(D).'*O.'*O))*inv(X.'*X);
   % grad = 2*X*D*inv(O)*inv(O).'*D.' - inv(X*D)*inv(X*D).'*O.'*O*inv(X*D)*inv(X*D).'*X*D*D.';
   grad = 2*X*D*inv(O)*inv(O).'*D.'-2*X*D*inv(O)*(inv(O).'*D.'*X.'*X*D*inv(O))^(-2)*inv(O).'*D.';
    div = zeros(size(X'))
    div(1,1) = 1e-4;
    compute_distortion(X'+div,Xbar')
    compute_distortion(X',Xbar')
end
