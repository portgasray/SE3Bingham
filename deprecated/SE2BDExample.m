%% Example for using the SE2 Bingham Distributioon

%% Parameters
C = [1, -0.1, 3 4;
     -0.1 1 5 6; 
      3 5  -1 0; 
      4 6  0 -1];

% C = [1, -0.1, 3 4;
%      -0.1 1 5 6; 
%       3 5  -1 0; 
%       4 6  0 -1];

%% Experiment with distribution
p = SE2BinghamDistribution(C);

% stocahstic sampling
s = p.sample(50000);
% s = p.sampleDeterministic

%% create plot of samples
rotateFirst = true;

% This algorithm assumes that the translation is performed after the
% rotation (Gilitschenski, 2014) 

% apply rotation to translation
% this is a different interpretation (first rotate, then translate)
phi = zeros(1, size(s,2));
translations = zeros(2, size(s,2));
if rotateFirst
    for i=1:size(s,2)
        phi(:,i) = 2 * atan2(s(2,i), s(1,i));
        R = [cos(phi(:,i)), -sin(phi(:,i)); sin(phi(:,i)), cos(phi(:,i))];
        translations(1,i) =  2 * (s(1,i) * s(3,i) - s(2,i)*s(4,i));
        translations(2,i) = 2 * (s(2,i) * s(3,i) + s(1,i)*s(4,i));
%         s(2:3, i) = R*s(2:3, i);
        translations(:,i) = R*translations(:,i);
    end
end

%% plot samples
figure(1)
quiver(translations(1,:), translations(2,:), cos(phi), sin(phi));
hold on
t_xy = zeros(2,1);
if rotateFirst
    m = p.mode();
    new_phi = 2 * atan2(m(2), m(1));
    R = [cos(new_phi), -sin(new_phi); sin(new_phi), cos(new_phi)];
    
    t_xy(1) =  2 * (m(1) * m(3) - m(2)*m(4));
    t_xy(2) = 2 * (m(2) * m(3) + m(1)*m(4));
    
    t_xy = R*t_xy;
    quiver(t_xy(1), t_xy(2), cos(new_phi),sin(new_phi), 'color','red', 'linewidth', 2);
% else
%     quiver(2 * (m(1) * m(3) - m(2)*m(4)), 2 * (m(2) * m(3) + m(1)*m(4)),cos(mu(1)),sin(mu(1)), 'color','red', 'linewidth', 2);
end
hold off
axis equal
xlabel('x')
ylabel('y')

%% plot 2D slice of density
% figure(2)
% step = 0.2;
% plotSize = 3;
% [x,y] = meshgrid([(mu(2)-plotSize):step:(mu(2)+plotSize)], [(mu(3)-plotSize):step:(mu(3)+plotSize)]);
% z = zeros(size(x));
% for i=1:size(x,1)
%     for j=1:size(x,2)
%         z(i,j) = p.pdf([mu(1), x(i,j), y(i,j)]');
%     end
% end
% surf(x,y,z)
% xlabel('x')
% ylabel('y')
% zlabel('f(\mu_1, x, y)');

