%% Example for using the SE3 Bingham Distributioon

%% Overall Settings
clear all;
addpath(genpath('E:\Repository\Experiments\Lei.Zhang\SE3Bingham\lib'));
addpath('E:\Repository\Experiments\Lei.Zhang\SE3Bingham\utils\export_fig');

N = 50;
fontSize = 14;
markerSize = 30;
lineWidth = 5;

%samples AutoScaleFactor
sampleSF = 0.5; % default 0.9
%mode AutoScaleFactor
modeSF = 1.2;

rotateFirst = true; %true %false
exportFig = false;
%% Parameters for the Bingham-Gaussian Distribution

C1 = diag([1 1 1 1]);

C11 = diag([1 1]);
C22 = [1 0.05; 0.05 1];
C33 = -diag([1 1]);
C2 = [C11 C22'; C22 C33];

C3 = -diag([1 1 1 1]);

C = [C1 C2'; C2 C3];


%% sample from SE3 Bingham distribution
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
theta1 = 2 * atan2(sqrt(s(2,:).^2 + s(3,:).^2 + s(4,:).^2), s(1,:));

% Rotation about an Arbitrary Vector
r_vec1 = [s(2,:)./sin(theta1/2); s(3,:)./sin(theta1/2); s(4,:)./sin(theta1/2);];
%% calculate the antipodal symmetry orientation -s(1:4,)
theta2 = 2 * atan2(sqrt(s(2,:).^2 + s(3,:).^2 + s(4,:).^2), -s(1,:));

r_vec2 = [-s(2,:)./sin(theta2/2); -s(3,:)./sin(theta2/2); -s(4,:)./sin(theta2/2)];

% A property: rad2deg(theta1 + theta2) = 360.0008, r_vec1 = -r_vec2  for each sample
%% calculate the translation from samples s(1:8,:)

% 2·(w0w5−w1w4+w2w7−w3w6)
% 2·(w0w6−w1w7−w2w4+w3w5)
% 2·(w0w7+w1w6−w2w5−w3w4)

% the sample s(1:4,:) from above
x = 2*(s(1,:).*s(6,:) - s(2,:).* s(5,:) + s(3,:).*s(8,:) - s(4,:).*s(7,:));
y = 2*(s(1,:).*s(7,:) - s(2,:).* s(8,:) + s(3,:).*s(5,:) - s(4,:).*s(6,:));
z = 2*(s(1,:).*s(8,:) - s(2,:).* s(7,:) + s(3,:).*s(6,:) - s(4,:).*s(5,:));
translation = [x; y; z];

% using the sample -s(1:4, :)
%% It is important to note, the translation will not change due to the Gaussian is assumed on dual part
%% solution 1: only nagete the first four -s(1:4,:)
% sym_t1 = 2*((-s(1,:)).*s(6,:) - (-s(2,:)).* s(5,:) + (-s(3,:)).*s(8,:) - (-s(4,:)).*s(7,:));
% sym_t2 = 2*((-s(1,:)).*s(7,:) - (-s(2,:)).* s(8,:) + (-s(3,:)).*s(5,:) - (-s(4,:)).*s(6,:));
% sym_t3 = 2*((-s(1,:)).*s(8,:) - (-s(2,:)).* s(7,:) + (-s(3,:)).*s(6,:) - (-s(4,:)).*s(5,:));
% sym_t = [sym_t1; sym_t2; sym_t3];

%% solution 2: nagete -s(1:8,:)
% sym_x = 2*((s(1,:)).*s(6,:) - (s(2,:)).* s(5,:) + (s(3,:).*s(8,:) - s(4,:).*s(7,:));
% sym_y = 2*((s(1,:)).*s(7,:) - (s(2,:)).* s(8,:) + (s(3,:)).*s(5,:) - s(4,:).*s(6,:));
% sym_z = 2*((s(1,:)).*s(8,:) - (s(2,:)).* s(7,:) + (s(3,:)).*s(6,:) - s(4,:).*s(5,:));
% sym_translation = [sym_x; sym_y; sym_z];

%% Rotaion first
% use cell to store the Rotation Matrices
Rotation = cell(size(s,2), 1);

if rotateFirst
    for i=1:size(s,2)
        Rotation{i} = angvec2r(theta1(i), r_vec1(:,i));
        translation(:,i) = Rotation{i}*translation(:,i) + [1;0;0];
    end
end

%% sample from s(1:8,:)

% samples generated from s(1:8, :)
% plot the rotational arbitrary vector

% offset setting for visualization, '0'-no offset
offset = 0;
figure(1);
% plot the postion by dots
scatter3(translation(1,:), translation(2,:), translation(3,:), markerSize, 'ob', 'fill');
hold on
% plot the orientaion by blue coloor
quiver3(translation(1,:), translation(2,:), translation(3,:), ...
    r_vec1(1,:), r_vec1(2,:), r_vec1(3,:), 'b', 'AutoScaleFactor', sampleSF, 'LineWidth', 1.2);
% hold on

% scatter3(translation(1,:)+offset, translation(2,:)+offset, translation(3,:)+offset, markerSize + 10, 'or'); %, 'fill'
% the same rotation represent by the r_vec2 with red color
% quiver3(translation(1,:)+offset, translation(2,:)+offset, translation(3,:)+offset, ...
%     r_vec2(1,:), r_vec2(2,:), r_vec2(3,:), 'r', 'AutoScaleFactor', sampleSF, 'LineWidth', 1.2);

%% convert the mode to angle and rotational arbitary vector and translation
% the second is the Gaussian distribution mode_gaussian = A_2 * M_bingham
% for -m is also the right answer ???

% rotation part: theta and vector - mode
m_theta1 = 2 * atan2(sqrt(m(2)^2 + m(3)^2 + m(4)^2), m(1));
m_rotaion_vec1 = [m(2)/sin(m_theta1/2), m(3)/sin(m_theta1/2), m(4)/sin(m_theta1/2)];

% translation part - mode
m_t1 = 2*(m(1)*m(6) - m(2)* m(5) + m(3)*m(8) - m(4)*m(7));
m_t2 = 2*(m(1)*m(7) - m(2)* m(8) + m(3)*m(5) - m(4)*m(6));
m_t3 = 2*(m(1)*m(8) - m(2)* m(7) + m(3)*m(6) - m(4)*m(5));
m_translation1 = [m_t1; m_t2; m_t3];

%% compute the antipodal one: -m( -m(1), )
a_m = -m;
% % rotation - antipodal mode
m_theta2 = 2 * atan2(sqrt(a_m(2)^2 + a_m(3)^2 + a_m(4)^2), a_m(1));
m_rotaion_vec2 = [a_m(2)/sin(m_theta2/2), a_m(3)/sin(m_theta2/2), a_m(4)/sin(m_theta2/2)];

% translation - antipodal mode
a_m_t1 = 2*(a_m(1)*a_m(6) - a_m(2)* a_m(5) + a_m(3)*a_m(8) - a_m(4)*a_m(7));
a_m_t2 = 2*(a_m(1)*a_m(7) - a_m(2)* a_m(8) + a_m(3)*a_m(5) - a_m(4)*a_m(6));
a_m_t3 = 2*(a_m(1)*a_m(8) - a_m(2)* a_m(7) + a_m(3)*a_m(6) - a_m(4)*a_m(5));

m_translation2 = [a_m_t1; a_m_t2; a_m_t3];
% important propety: m_translation2 =  m_translation1
%% plot modes
% logic: if rotation happends, 
hold on
if rotateFirst
    R_mode = angvec2r(m_theta1, m_rotaion_vec1);
    m_tNew = R_mode*m_translation1;
    quiver3(m_tNew(1),m_tNew(2),m_tNew(3), ...
        m_rotaion_vec1(1), m_rotaion_vec1(2), m_rotaion_vec1(3), ...
        'r', 'LineWidth', lineWidth, 'AutoScaleFactor',1.2);
    scatter3(m_tNew(1), m_tNew(2), m_tNew(3), markerSize+10, 'or', 'fill');
%     R_mode_anti = angvec2r(m_theta2, m_rotaion_vec2);
%     rotated_mode_anti = R_mode_anti*m_translation2;
%     hold on;
%     quiver3(rotated_mode_anti(1),rotated_mode_anti(2),rotated_mode_anti(3), ...
%         m_rotaion_vec2(1), m_rotaion_vec2(2), m_rotaion_vec2(3), ...
%         'b', 'LineWidth', lineWidth, 'AutoScaleFactor',1.2);
else
    scatter3(m_translation1(1), m_translation1(2), m_translation1(3), markerSize+10, 'ob', 'fill');
    quiver3(m_translation1(1), m_translation1(2), ...
        m_translation1(3), m_rotaion_vec1(1), m_rotaion_vec1(2), m_rotaion_vec1(3), ...
        'b', 'LineWidth', lineWidth, 'AutoScaleFactor', modeSF);
    hold on;
    quiver3(m_translation2(1), m_translation2(2), m_translation2(3), ...
        m_rotaion_vec2(1), m_rotaion_vec2(2), m_rotaion_vec2(3), ...
        'r', 'LineWidth', lineWidth, 'AutoScaleFactor', modeSF);
end
hold off;




%% second plot for rotated from the same samples
% when rotaion first; use cell to store the Rotation Matrices
% Rotation = cell(size(s,2), 1);
% rotateFirst = true; %true %false
% if rotateFirst
%     for i=1:size(s,2)
%         Rotation{i} = angvec2r(theta1(i), v(i,:));
%         t(:,i) = Rotation{i}*t(:,i);
%     end
% end
% 
% offset = 0;
% figure(2);
% scatter3(t(1,:), t(2,:), t(3,:), markerSize-20, 'ob', 'fill');
% hold on
% scatter3(t1+offset, t2+offset, t3+offset, markerSize, 'or'); %, 'fill'

% quiver3(t(1,:), t(2,:), t(3,:), v1, v2, v3, 'b', 'AutoScaleFactor', 0.3);
% 
% hold on
% 
% quiver3(t1+offset, t2+offset, t3+offset, sym_v1, sym_v2, sym_v3, 'r', 'AutoScaleFactor', 0.5);

% hold off
%% configuration for plot
set(gca, 'FontSize', fontSize);
axis vis3d; axis equal;

% axes.SortMethod='ChildOrder';
xlabel('x'); ylabel('y'); zlabel('z');
grid on; grid minor;
box on;

% az = -47;
% el = 9;
% view(az, el);

% set(findall(gcf,'type','text'),'FontSize',fontSize);

filename = ['se3bingham-pose-' num2str(N) '.pdf'];
% filename = sprintf('se3bingham-pose.pdf');

if rotateFirst
%     filename = sprintf('se3bingham-pose-transformed.pdf');
    filename = ['se3bingham-pose-' num2str(N) '-transformed.pdf'];
end

%% export figure

if exportFig
%     filename = sprintf('se3bingham-pose-transformed.pdf');
    export_fig(gcf, filename, '-transparent');
end
