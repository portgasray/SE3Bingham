function plot_3d_projection()
%PLOT_3D_PROJECTION Summary of this function goes here
%   for R^{4}: using 4 S^{3} to display

[SX,SY,SZ] = sphere(150); % Create smooth sphere
n = size(SX,1);

C1 = zeros(n);
C2 = zeros(n);
C3 = zeros(n);
C4 = zeros(n);

%BinghamDistribution.pdf(this, xa)
% Parameters:
%   xa (d x n matrix)
%       each column represents one of the n points in R^d that the
%       pdf is evaluated at; can be just a (d x 1) vector as well
% Returns:
%   p (1 x n row vector)
%       values of the pdf at each column of xa

% test pdf
%rng default
%   testpoints = rand(B.dim, 20);
%   testpoints = bsxfun(@rdivide,testpoints, sum(testpoints));
%   for j=1:size(testpoints,2)
%       testCase.verifyEqual(B.pdf(testpoints(:,j)), 1/B.F*exp(testpoints(:,j)'*M*diag(Z)*M'*testpoints(:,j)), 'RelTol', 1E-10);
%   end

for i=1:n
   for j=1:n
      q1 = [0 SX(i,j) SY(i,j) SZ(i,j)];
      q2 = [SX(i,j) 0 SY(i,j) SZ(i,j)];
      q3 = [SX(i,j) SY(i,j) 0 SZ(i,j)];
      q4 = [SX(i,j) SY(i,j) SZ(i,j) 0];
%       C1(i,j) = BinghamDistribution.pdf(q1, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
%       C2(i,j) = BinghamDistribution.pdf(q2, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
%       C3(i,j) = BinghamDistribution.pdf(q3, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
%       C4(i,j) = BinghamDistribution.pdf(q4, Z(1), Z(2), Z(3), V(:,1), V(:,2), V(:,3), F);
      C
   end
end

C1 = C1./max(max(C1));
C2 = C2./max(max(C2));
C3 = C3./max(max(C3));
C4 = C4./max(max(C4));

subplot(2,2,1);
surf(SX,SY,SZ,C1, 'EdgeColor', 'none');
axis vis3d;
subplot(2,2,2);
surf(SX,SY,SZ,C2, 'EdgeColor', 'none');
axis vis3d;
subplot(2,2,3);
surf(SX,SY,SZ,C3, 'EdgeColor', 'none');
axis vis3d;
subplot(2,2,4);
surf(SX,SY,SZ,C4, 'EdgeColor', 'none');
axis vis3d;

end

