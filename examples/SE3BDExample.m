%% Example for using the SE3 Bingham Distributioon

%% Overall Settings
clear all;
addpath(genpath('../lib'));
addpath('../utils/export_fig');

fontSize = 14;
markerSize = 10;
lineWidth = 5;

rotateFirst = false; %true %false
%% Parameters for the Bingham-Gaussian Distribution
% group 1
C1 = -diag([2 2 2 2]);

C11 = -diag([2 3]);
C22 = [0.1 0.2; 0.01 0.3];
C33 = -diag([1 2]);
C2 = [C11 C22'; C22 C33];

C3 = -diag([2 2 2 2]);

C = [C1 C2'; C2 C3];

% group 2
% C1 = -diag([1 1 1 1]);
% 
% C11 = -diag([1 1]);
% C22 = [0.1 0.1; 0.1 0.1];
% C33 = -diag([0.1 0.1]);
% C2 = [C11 C22'; C22 C33];
% 
% C3 = -diag([1 1 1 1]);
% 
% C = [C1 C2'; C2 C3];

%% calculating the mode (deprecated)
% binghamC = C1 - C2'*pinv(C3)*C2;
% [M, Z] = eig(binghamC);
% 
% [Z,order] = sort(diag(Z),'ascend');
% M = M(:,order);
% 
% m = zeros(8,1);
% m(1:4) = M(:,4);
% m(5:8) = - pinv(C3)*C2*m(1:4);

%% sample from SE3 Bingham distribution
% [M, Z] = eig(C2);
% Z = sort(diag(Z), 'ascend');
% Z(end) = 0;
% b = BinghamDistribution([-30 -10 -3 0]', eye(4,4));
% s = b.sample(N);

N = 10;
p = SE3BinghamDistribution(C);
s = p.sample(N);
% mode 8-length vector
m = p.mode();

%%  Algorithm 1 recover the rotation angle, arbitray vector and translation from samples s
% This algorithm assumes that the translation is performed after the
% rotation (apply rotation to translation  )
% this is a different interpretation (first rotate, then translate)

% s(1:4) represent for rotation thus s(1:4) and -s(1:4) represent the same
% orientation

% s(5:8) is the dual part, which can be converted to translation in 3D

%% calculate the rotation from sample s(1:4,:)
theta = 2 * atan2(sqrt(s(2,:).^2 + s(3,:).^2 + s(4,:).^2), s(1,:));

v1 = s(2,:)./sin(theta/2);
v2 = s(3,:)./sin(theta/2);
v3 = s(4,:)./sin(theta/2);
v = [v1; v2; v3]';

%% calculate the antipodal symmetry orientation -s(1:4,)
% sym_theta = 2 * atan2(sqrt(s(2,:).^2 + s(3,:).^2 + s(4,:).^2), -s(1,:));
% 
% sym_v1 = -s(2,:)./sin(sym_theta/2);
% sym_v2 = -s(3,:)./sin(sym_theta/2);
% sym_v3 = -s(4,:)./sin(sym_theta/2);
% sym_v = [sym_v1; sym_v2; sym_v3]';

% nice properties: rad2deg(theta + sym_theta) = 360.0008, sym_v = -v  for each sample
%% calculate the translation from samples s(1:8,:)

% 2·(w0w5−w1w4+w2w7−w3w6)
% 2·(w0w6−w1w7−w2w4+w3w5)
% 2·(w0w7+w1w6−w2w5−w3w4)

% the sample s(1:4, :) is from above
t1 = 2*(s(1,:).*s(6,:) - s(2,:).* s(5,:) + s(3,:).*s(8,:) - s(4,:).*s(7,:));
t2 = 2*(s(1,:).*s(7,:) - s(2,:).* s(8,:) + s(3,:).*s(5,:) - s(4,:).*s(6,:));
t3 = 2*(s(1,:).*s(8,:) - s(2,:).* s(7,:) + s(3,:).*s(6,:) - s(4,:).*s(5,:));
t = [t1; t2; t3];

% using the sample -s(1:4, :)
% sym_t1 = 2*((-s(1,:)).*s(6,:) - (-s(2,:)).* s(5,:) + (-s(3,:)).*s(8,:) - (-s(4,:)).*s(7,:));
% sym_t2 = 2*((-s(1,:)).*s(7,:) - (-s(2,:)).* s(8,:) + (-s(3,:)).*s(5,:) - (-s(4,:)).*s(6,:));
% sym_t3 = 2*((-s(1,:)).*s(8,:) - (-s(2,:)).* s(7,:) + (-s(3,:)).*s(6,:) - (-s(4,:)).*s(5,:));
% sym_t = [sym_t1; sym_t2; sym_t3];

% solution 2: nagete -s(1:8,:)
% sym_t1 = 2*((s(1,:)).*s(6,:) - (s(2,:)).* s(5,:) + (s(3,:).*s(8,:) - s(4,:).*s(7,:));
% sym_t2 = 2*((s(1,:)).*s(7,:) - (s(2,:)).* s(8,:) + (s(3,:)).*s(5,:) - s(4,:).*s(6,:));
% sym_t3 = 2*((s(1,:)).*s(8,:) - (s(2,:)).* s(7,:) + (s(3,:)).*s(6,:) - s(4,:).*s(5,:));
% sym_t = [sym_t1; sym_t2; sym_t3];


%% when rotaion first (deprecated)
% if rotateFirst
%     for i=1:size(s,2)
%         R = angvec2r(theta(i), v(i,:));
%         t(:,i) = R*t(:,i);
%     end
% end

%% when rotaion first; use cell to store the Rotation Matrices
Rotation = cell(size(s,2), 1);
if rotateFirst
    for i=1:size(s,2)
        Rotation{i} = angvec2r(theta(i), v(i,:));
        t(:,i) = Rotation{i}*t(:,i);
    end
end

%% plot samples
clf
hold on
% samples generated from s(1:8, :)
% plot the rotational arbitrary vector
% quiver3(t1, t2, t3, v1, v2, v3, 'b', 'AutoScaleFactor', 0.3);
% plot the postion by dots
scatter3(t1, t2, t3, markerSize, 'ob');


% quiver3(sym_t1, sym_t2, sym_t3, sym_v1, sym_v2, sym_v3, 'r', 'AutoScaleFactor', 0.3);
% scatter3(sym_t1, sym_t2, sym_t3, markerSize, 'or', 'fill');

%% plot modes

%% convert the mode to angle and rotational arbitary vector and translation
% the second is the Gaussian distribution mode_gaussian = A_2 * M_bingham
% for -m is also the right answer

% rotation part: theta and vector - mode
m_theta = 2 * atan2(sqrt(m(2)^2 + m(3)^2 + m(4)^2), m(1));
m_vec = [m(2)/sin(m_theta/2), m(3)/sin(m_theta/2), m(4)/sin(m_theta/2)];

% translation part - mode
m_t1 = 2*(m(1)*m(6) - m(2)* m(5) + m(3)*m(8) - m(4)*m(7));
m_t2 = 2*(m(1)*m(7) - m(2)* m(8) + m(3)*m(5) - m(4)*m(6));
m_t3 = 2*(m(1)*m(8) - m(2)* m(7) + m(3)*m(6) - m(4)*m(5));
m_t = [m_t1; m_t2; m_t3];

% compute the antipodal one: -m( -m(1), )
% a_m = -m;
% 
% 
% % rotation - antipodal mode
% a_m_theta = 2 * atan2(sqrt(a_m(2)^2 + a_m(3)^2 + a_m(4)^2), a_m(1));
% a_m_vec = [a_m(2)/sin(a_m_theta/2), a_m(3)/sin(a_m_theta/2), a_m(4)/sin(a_m_theta/2)];
% 
% % translation - antipodal mode
% a_m_t1 = 2*(a_m(1)*a_m(6) - a_m(2)* a_m(5) + a_m(3)*a_m(8) - a_m(4)*a_m(7));
% a_m_t2 = 2*(a_m(1)*a_m(7) - a_m(2)* a_m(8) + a_m(3)*a_m(5) - a_m(4)*a_m(6));
% a_m_t3 = 2*(a_m(1)*a_m(8) - a_m(2)* a_m(7) + a_m(3)*a_m(6) - a_m(4)*a_m(5));
% 
% a_m_t = [a_m_t1; a_m_t2; a_m_t3];

% hold on
% if rotateFirst
%     R_mode = angvec2r(mode_theta, mode_v);
%     m_tNew = R_mode*m_t;
%     quiver3(m_tNew(1),m_tNew(2),m_tNew(3), mode_v(1), mode_v(2), mode_v(3), 'c', 'LineWidth', lineWidth);
% else
%     quiver3(m_t(1), m_t(2), m_t(3), mode_v(1), mode_v(2), mode_v(3), 'b', 'LineWidth', lineWidth);
%     hold on
%     quiver3(anti_m_t(1), anti_m_t(2), anti_m_t(3), anti_m_vec(1), anti_m_vec(2), anti_m_vec(3), 'r', 'LineWidth', lineWidth);
% end
% hold off

%% configuration for plot
set(gca, 'FontSize', fontSize);
axis vis3d; axis equal;

% axes.SortMethod='ChildOrder';
xlabel('x'); ylabel('y'); zlabel('z');
grid on; grid minor;
box on;
% view(3);
set(findall(gcf,'type','text'),'FontSize',fontSize);


filename = sprintf('se3-bingham-plot3-raw.pdf');
if rotateFirst
    filename = sprintf('se3-bingham-plot3-rotated.pdf');
end

%% export figure
% export_fig(gcf, filename, '-transparent');

%% plot 3D slice of density (de)
% figure(2)
% step = 0.2;
% plotSize = 3;
% [x,y] = meshgrid([(m(7)-plotSize):step:(m(3)+plotSize)], [(m(8)-plotSize):step:(m(4)+plotSize)]);
% z = zeros(size(x));
% for i=1:size(x,1)
%     for j=1:size(x,2)
% %         z(i,j) = p.pdf([m(5), m(6), x(i,j), y(i,j)]');
%     end
% end
% surf(x,y,z)
% xlabel('x')
% ylabel('y')
% zlabel('f(\mu_1, x, y)');
% 
% set(gca, 'FontSize', fontSize);
% axis vis3d; axes.SortMethod='ChildOrder';
% grid on; grid minor;
% 
% filename = sprintf('se3-bingham-sample-5000-2d-density.pdf');
% if ~rotateFirst
% %     export_fig(gcf, filename, '-transparent');
% end

