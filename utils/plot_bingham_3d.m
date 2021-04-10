% b1 = BinghamDistribution([-1 -1 0]', eye(3,3));
% b2 = BinghamDistribution([-5 -1 0]', eye(3,3));
% b3 = BinghamDistribution([-50 -1 0]', eye(3,3));
% figure(1);
% b1.plot;
% figure(2);
% b2.plot;
% figure(3);
% b3.plot;



% [X,Y,Z] = peaks(25);
% colormap winter;
% CO(:,:,1) = zeros(25).*linspace(1,0,25); % red
% CO(:,:,2) = ones(25).*linspace(0,0,25); % green
% CO(:,:,3) = ones(25).*linspace(1,0,25); % blue
% surf(X,Y,Z,CO);
% % surf(X,Y,Z, [1 1 0], 'none','EdgeColor','none')
% colorbar;

% 
% [X,Y] = meshgrid(1:0.5:10,1:20);
% Z = sin(X) + cos(Y);
% C = '';
% surf(X,Y,Z,C)
% colorbar

% tiledlayout(1,2)
% ax1 = nexttile;
% surf(peaks);
% shading interp;
% colormap(ax1,parula(10));

ax2 = nexttile;
surf(peaks);
shading interp;
colormap(ax2,cool(5));