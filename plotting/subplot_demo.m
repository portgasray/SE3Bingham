% subplot(2,2,1)
% x = linspace(0,10);
% y1 = sin(x);
% plot(x,y1)
% title('Subplot 1: sin(x)')
% 
% subplot(2,2,2)
% y2 = sin(2*x);
% plot(x,y2)
% title('Subplot 2: sin(2x)')
% 
% subplot(2,2,3)
% y3 = sin(4*x);
% plot(x,y3)
% title('Subplot 3: sin(4x)')
% 
% subplot(2,2,4)
% y4 = sin(8*x);
% plot(x,y4)
% title('Subplot 4: sin(8x)');

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


% if nargin<3; gap = .02; end
% if nargin<4 || isempty(marg_h); marg_h = .05; end
% if nargin<5; marg_w = .05; end
% 
% if numel(gap)==1; 
%     gap = [gap gap];
% end
% if numel(marg_w)==1; 
%     marg_w = [marg_w marg_w];
% end
% if numel(marg_h)==1; 
%     marg_h = [marg_h marg_h];
% end
% 
% axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
% axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
% 
% py = 1-marg_h(2)-axh; 
% 
% ha = zeros(Nh*Nw,1);
% ii = 0;
% for ih = 1:Nh
%     px = marg_w(1);
% 
%     for ix = 1:Nw
%         ii = ii+1;
%         ha(ii) = axes('Units','normalized', ...
%             'Position',[px py axw axh], ...
%             'XTickLabel','', ...
%             'YTickLabel','');
%         px = px+axw+gap(2);
%     end
%     py = py-axh-gap(1);
% end


fontSize = 12;
addpath('E:\Repository\Experiments\Lei.Zhang\utils\export_fig');
% [X,Y] = meshgrid(-5:.5:5);
% Z = Y.*sin(X) - X.*cos(Y);
% s = surf(X,Y,Z,'FaceAlpha',0.5)
% set(gca, 'FontSize', fontSize);
% axis vis3d; axis equal;
% xlabel('x');ylabel('y');zlabel('z');
% grid on;  grid minor;
% box on;
% % filename = ['projected_bingham_' num2str(i) '.pdf'];
% filename = sprintf('demo_for+faceAlpha.pdf');
% export_fig(gcf, filename, '-transparent');

colormap(cool)
surf(peaks)
% export_fig('output.eps','-eps','-pdf','-opengl','-m3','-q101','-transparent','-nocrop');

%generate the EPS export using the following command:
print(gcf, '-painters', '-loose', '-depsc2', 'test.eps')