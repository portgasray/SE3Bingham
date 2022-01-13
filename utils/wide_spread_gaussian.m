
addpath('E:\Repository\Experiments\Lei.Zhang\SE3Bingham\utils\export_fig\');

%% 1D gaussian
x = [-3:.1:3];
y = normpdf(x,0,1);
plot(x,y)
%% wide spread 3D gaussian distribution
% mu = [0 0];
% Sigma = [0.9 0.2; 0.2 0.9];
% 
% x1 = -3:0.2:3;
% x2 = -3:0.2:3;
% [X1,X2] = meshgrid(x1,x2);
% X = [X1(:) X2(:)];
% 
% y = mvnpdf(X,mu,Sigma);
% y = reshape(y,length(x2),length(x1));

%% surf plot 3d
% surf(x1,x2,y)
% caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
% axis([-3 3 -3 3 0 0.4])
% xlabel('x1')
% ylabel('x2')
% zlabel('Probability Density')

%% export to pdf
set(gca, 'FontSize', 16);

% set(gca,'visible','off')
% axis vis3d; axis equal;
% xlabel('x');ylabel('y');zlabel('z');
axis off
grid on; 
grid minor;
box on;


% set(findall(gcf,'type','text'),'FontSize',16);
filename = ['wide_spread_gaussian.pdf'];
export_fig(gcf, filename, '-transparent');
