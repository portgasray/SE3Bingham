clear all;
% syms d0 d1 d2 d3 d4 d5 d6 d7 real
% r = [d0 d1 d2 d3];
% d = [d4 d5 d6 d7];

% 2D planar
syms d0 d1 d2 d3 real
xs = [d0 d1];
xt = [d2 d3];

C1 = -diag([2 3]);
C2 = [0.1 0.2; 0.01 0.3];
C3 = -diag([1 2]);
            
C = [C1 C2'; C2 C3];


% seA = SE2BinghamDistribution(C);
% seB = SE2BinghamDistribution(C1,C2,C3);

T1 = C1-(C2' * pinv(C3)*C2);
T2 = - pinv(C3)*C2;
new_x = xt' - (T2*xs');
norm = norm(new_x);
norm

%read file - S3 data
% syn_cube_all in riss-bingham
X = read_matrix('syn_cube_all.txt'); % file under riss-bingham matlab folder
X = X'; % change size to (4, n_samples)
B = (X');
% plot_bingham_3d_projections for 4D-unit vector S3

bingham_fit_ransac(X);