% % syn_cube_all in riss-bingham
% % X = read_matrix('test_sampled_quat_trimed.txt'); % file under riss-bingham matlab folder
% X = read_matrix('syn_cube_all.txt');
% 
% dim = size(X);
% N = dim(1);
% d = dim(2);
% angles = zeros(N, 1);
% translations = zeros(N,2);
% for i=1:N
%     angles(i) = 2 * atan2(X(i,2), X(i,1));
%     translations(i,1) =  2 * (X(i,1) * X(i,3) - X(i,2)*X(i,4));
%     translations(i,2) = 2 * (X(i,2) * X(i,3) + X(i,1)*X(i,4));
% end
% figure(1)
% % % quiver(s(2,:),s(3,:),cos(s(1,:)),sin(s(1,:)));
% quiver(translations(:,1), translations(:,2),cos(angles),sin(angles));


%% March, 28, 2021
X1 =  readmatrix('E:\Repository\Experiments\Lei.Zhang\data\test_sampled_quat.txt');
quat = X1(:,1:4)';
q = quat(:,1);

M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
Z = [-10 -2 -1 0]';
B = BinghamDistribution(Z,M);
N = 1000;
% samples = B_quat.sample(N);
samples = B_quat.sampleGlover(N);
X = samples';

dim = size(X);
N = dim(1);
d = dim(2);
angles = zeros(N, 1);
translations = zeros(N,2);
for i=1:N
    angles(i) = 2 * atan2(X(i,2), X(i,1));
    translations(i,1) =  2 * (X(i,1) * X(i,3) - X(i,2)*X(i,4));
    translations(i,2) = 2 * (X(i,2) * X(i,3) + X(i,1)*X(i,4));
end
figure(1)
% % quiver(s(2,:),s(3,:),cos(s(1,:)),sin(s(1,:)));
quiver(translations(:,1), translations(:,2),cos(angles),sin(angles));