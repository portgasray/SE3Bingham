%% Example for using the SE2 Bingham Distributioon
addpath(genpath('E:\Repository\Experiments\Lei.Zhang\SE3Bingham\lib'));
%% Parameters

% group 1
% C1 = -diag([2 3]);
% C2 = [0.1 0.2; 0.01 0.3];
% C3 = -diag([1 2]);

C1 = -diag([2 3]);
C2 = [0.1 0.2; 0.01 3];
C3 = -diag([1 0.2]);

C = [C1 C2'; C2 C3];

%% marginal distributio of Bingham Distribution
binghamC = C1 - C2'*pinv(C3)*C2;
[M, Z] = eig(binghamC);

% [Z,order] = sort(diag(Z),'ascend');
% M = M(:,order);

if Z(1,1)>Z(2,2)
    M = [M(:, 2) M(:,1)];
    Z = diag([Z(2,2) Z(1,1)]);
end

N = 5000;
%% Experiment with distribution
p = SE2BinghamDistribution(C);
s = p.sample(N);
m = p.mode();
%% sample from Bingham distribution
% [M, Z] = eig(C2);
% Z = sort(diag(Z), 'ascend');
% Z(end) = 0;
% b = BinghamDistribution([-30 -10 -3 0]', eye(4,4));
% s = b.sample(N);
%% create plot of samples
rotateFirst = true; %true %false

% This algorithm assumes that the translation is performed after the
% rotation (Gilitschenski, 2014) 

% apply rotation to translation
% this is a different interpretation (first rotate, then translate)
theta = 2 * atan2(s(2,:), s(1,:));

if rotateFirst
    for i=1:size(s,2)
        R = [cos(theta(i)), -sin(theta(i)); sin(theta(i)), cos(theta(i))];
        s(3,i) =  2 * (s(1,i) * s(3,i) - s(2,i)*s(4,i));
        s(4,i) = 2 * (s(2,i) * s(3,i) + s(1,i)*s(4,i));
        s(3:4, i) = R*s(3:4, i);
    end
end

%% plot samples
figure(1)
quiver(s(3,:),s(4,:), cos(theta), sin(theta));

%% plot sample mode
hold on
mTheta =  2 * atan2(m(2), m(1));
if rotateFirst
    R = [cos(mTheta), -sin(mTheta); sin(mTheta), cos(mTheta)];
    mNew(3:4) = R*m(3:4);
    quiver(mNew(3),mNew(4),cos(mTheta),sin(mTheta), 'color','red', 'linewidth', 2);
else
    quiver(m(3),m(4),cos(mTheta),sin(mTheta), 'color','red', 'linewidth', 2);
    quiver(-m(3),-m(4),cos(mTheta),sin(mTheta), 'color','red', 'linewidth', 2);
end
hold off

%% configuration for plot to export
fontSize = 12;
set(gca, 'FontSize', fontSize);
% axis vis3d; 
axis equal; axes.SortMethod='ChildOrder';
xlabel('x'); ylabel('y'); zlabel('z');
grid on; grid minor;
box on;

set(findall(gcf,'type','text'),'FontSize',fontSize);
addpath('E:\Repository\Experiments\Lei.Zhang\utils\export_fig');

filename = sprintf('se2-bingham-sample-5000.pdf');
if rotateFirst
    filename = sprintf('se2-bingham-sample-5000-rotated.pdf');
end
% export_fig(gcf, filename, '-transparent');

%% plot 2D slice of density
figure(2)
step = 0.2;
plotSize = 3;
[x,y] = meshgrid([(m(3)-plotSize):step:(m(3)+plotSize)], [(m(4)-plotSize):step:(m(4)+plotSize)]);
z = zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        z(i,j) = p.pdf([m(1), m(2), x(i,j), y(i,j)]');
    end
end
surf(x,y,z)
xlabel('x')
ylabel('y')
zlabel('f(\mu_1, x, y)');

set(gca, 'FontSize', fontSize);
axis vis3d; axes.SortMethod='ChildOrder';
grid on; grid minor;

filename = sprintf('se2-bingham-sample-5000-2d-density.pdf');
if ~rotateFirst
%     export_fig(gcf, filename, '-transparent');
end
