addpath(genpath('E:\Workspace\MATLAB\mtex-5.5.0'));

%% Parameters for SE3Bingham Dist
C1 = -diag([2 2 2 2]);

C11 = -diag([2 3]);
C22 = [0.1 0.2; 0.01 0.3];
C33 = -diag([0.1 2]);

C2 = [C11 C22'; C22 C33];

C3 = -diag([2 2 2 2]);

% C2 = [0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1];
C = [C1 C2'; C2 C3];

binghamC = C1 - C2'*pinv(C3)*C2;
[M, Z] = eig(binghamC);

[Z,order] = sort(diag(Z),'ascend');
M = M(:,order);

z=diag(Z)-Z(2,2);

%% odf plot
cs = crystalSymmetry('1');

% kappa = [100 90 80 0];   % shape parameters
% U     = eye(4);          % orthogonal matrix

odf = BinghamODF(z,M,cs);

h = [Miller(0,0,1,cs) Miller(1,0,0,cs) Miller(1,1,1,cs)];
plotPDF(odf,h,'antipodal','silent');


%% parameter form Glover's PhD.Thesis
% kappa = [-43.2, -2.3 0 0];

% kappa = [-49.0, -1.0 0 0];

% kappa = [-100, 0 0 0];
% 
% kappa = [-29.9, -4.7 0 0];

% kappa = [-41.6, -3.3 0 0];

% kappa = [-94.0, -0.6 0 0];

% kappa = [-73.6 -1.8 0 0];

kappa = [-32.8 -3.7 0 0];

U = eye(4);
odf = BinghamODF(kappa,U,cs);
h = [Miller(0,0,1,cs) Miller(1,0,0,cs) Miller(1,1,1,cs)];
plotPDF(odf,h,'antipodal','silent');