clear all;
addpath(genpath('E:\Workspace\MATLAB\bingham'));

addpath('E:\Repository\Experiments\Lei.Zhang\data');
addpath(genpath('E:\Workspace\MATLAB\riss_bingham\'));

addpath('E:\Repository\Experiments\Lei.Zhang\SE3Bingham\utils\export_fig\');

exportFig = true;


%%first input
% X = read_matrix('trimed_test_sampled_quat.txt');
% quat = quaternion(X);
X = readmatrix('syn_cube_all.txt');
B = bingham_fit(X);

%%second input
% B.d = 4;
% %diagonal matrix
% % B.Z = [-30,-30,-600];
% 
% % B.Z = [-60,-60,-900];
% 
% B.Z = [-600,-600,-900];
% %orthogonal matrix
% B.V = eye(4);
% [B.F B.df] = bingham_F(B.Z);


%% third input 
C1 = -diag([2 2 2 2]);

C11 = -diag([2 3]);
C22 = [0.1 0.2; 0.01 0.3];
C33 = -diag([0.1 2]);
C2 = [C11 C22'; C22 C33];

C3 = -diag([2 2 2 2]);

C = [C1 C2'; C2 C3];

binghamC = C1 - C2'*pinv(C3)*C2;
[M, Z] = eig(binghamC);

[Z,order] = sort(diag(Z),'ascend');
M = M(:,order);
z=Z-Z(end);

% B.d = 4;
% B.Z = z(1:3)';
% [B.F B.df] = bingham_F(B.Z);
%%
V = B.V; Z = B.Z; F = B.F;

%%
% implement the content in plot_bingham_3d_projections(B.V, B.Z, B.F);
%clf;
[SX,SY,SZ] = sphere(30);
n = size(SX,1);

C1 = zeros(n);
C2 = zeros(n);
C3 = zeros(n);
C4 = zeros(n);

for i=1:n
   for j=1:n
      q1 = [0 SX(i,j) SY(i,j) SZ(i,j)];
      q2 = [SX(i,j) 0 SY(i,j) SZ(i,j)];
      q3 = [SX(i,j) SY(i,j) 0 SZ(i,j)];
      q4 = [SX(i,j) SY(i,j) SZ(i,j) 0];
%       To compute the PDF of a unit vector x under a Bingham B, use `f = bingham_pdf(x,B);`
% Here, I computer 4 
      C1(i,j) = bingham_pdf_3d(q1, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
      C2(i,j) = bingham_pdf_3d(q2, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
      C3(i,j) = bingham_pdf_3d(q3, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
      C4(i,j) = bingham_pdf_3d(q4, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
   end
end

C1 = C1./max(max(C1));
C2 = C2./max(max(C2));
C3 = C3./max(max(C3));
C4 = C4./max(max(C4));

fontSize = 16;
% N = 2;
% p = panel();
% 
% p.select('all');



%% plot
% surf(SX,SY,SZ,'EdgeColor', 'none','FaceAlpha', .3);
colormap(.5*cool+.5);

% subplot(2,2,1);
% surf(SX,SY,SZ,C1, 'EdgeColor', 'none','FaceAlpha', .5);
% axis vis3d; axis equal;
% xlabel('x');ylabel('y');zlabel('z');
% grid on; grid minor;
% % box on;
% 
% p(1, 1).select();
% subplot(2,2,2);
% % surf(SX,SY,SZ,C2, 'EdgeColor', 'none','FaceAlpha', .3);
% surf(SX,SY,SZ,C2, 'EdgeColor', 'none','FaceAlpha', .3);
% axis vis3d; axis equal;  % axis tight; axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% grid on; grid minor;
% 
% subplot(2,2,3);
% surf(SX,SY,SZ,C3, 'EdgeColor', 'none','FaceAlpha', .3);
% axis vis3d; axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% grid on; grid minor;
% 
% subplot(2,2,4);
% surf(SX,SY,SZ,C4, 'EdgeColor', 'none','FaceAlpha', .3);
% axis vis3d; axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% grid on; grid minor;

% colorbar;

%% cell to store matrices
%  n = 10 ;
%  M = cell(n, 1) ;
%  for k = 1 : n
%     M{k} = 20*k + rand(10) ;
%  end

M = cell(4, 1);
M{1} = C1;
M{2} = C2;
M{3} = C3;
M{4} = C4;

for i=1:4
    figure(i);
%     p = panel();
%     p.pack('h', {95 []});
% %     
%     h_axis = p(1).select();
%     s = surf(SX, SY, SZ, M{i}, 'EdgeColor', 'none', 'FaceAlpha', .5); 
    s = surf(SX, SY, SZ, M{i}, 'EdgeColor', 'none', 'FaceAlpha', .7); 
    %'interp'
%     
%     % sometimes you'll want to use some other function than
%     % Panel to create one or more axes. for instance,
%     % colorbar...
%     h_colorbar_axis = colorbar('peer', h_axis);
%     p(2).select(h_colorbar_axis);
%     axis([0 100 -3 3]);

%     colorbar;
%     colormap(.5*cool+.5);
    
%     s.EdgeColor='flat' %interp
%     s.FaceColor = 'interp';
%     s.LineStyle='none';
%     if i == 1;
%         surf(SX,SY,SZ, C1, 'EdgeColor', 'none');
%     end
%     if i == 2;
%         surf(SX,SY,SZ, C2, 'EdgeColor', 'none');
%     end
%     if i == 3;
%         surf(SX,SY,SZ, C3, 'EdgeColor', 'none');
%     end
%     if i == 4;
%         surf(SX,SY,SZ, C4, 'EdgeColor', 'none');
%     end
    % configuration for plot
    set(gca, 'FontSize', fontSize);
    axis vis3d; 
    axis equal;
%     axis off;
    xlabel('x');ylabel('y');zlabel('z');
    grid on; 
    grid minor;
    box on;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(findall(gcf,'type','text'),'FontSize',fontSize);
    filename = ['projected_bingham_' num2str(i) '.pdf'];
   
    if exportFig;
        export_fig(gcf, filename, '-transparent');
    end
    
%     export_fig filename -transparent -painters
end

% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% filename = sprintf('3-sphere-projected_bingham.pdf');
% set(findall(gcf,'type','text'),'FontSize',fontSize);
% export_fig(gcf, filename, '-transparent');

