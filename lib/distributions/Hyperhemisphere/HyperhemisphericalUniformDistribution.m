classdef HyperhemisphericalUniformDistribution < AbstractHyperhemisphericalDistribution & AbstractUniformDistribution
    methods
        function this = HyperhemisphericalUniformDistribution(dim_)
            assert(dim_>=2)
            this.dim = dim_;
        end
        
        function vals = pdf(this, xa)
            % Evaluates pdf at each column of xa
            % Parameters:
            %   xa (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of xa
            
            % pdf is always 1 divided by surface area (integrates to 1).
            vals = (2/AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim))*ones(1,size(xa,2));
        end
        
        function s = sample(this, n)
            % Stocahastic sampling
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column     
            %
            % General algorithm based on "A note on a method for generating 
            % points uniformly on n-dimensional spheres" by Mervin E. Muller April 1959.
            % Algorithm for 2-sphere based on "Spherical sampling by archimedes' theorem"
            % by Min-Zhi Shao and Norman Badler, 1996

            assert(isscalar(n));
            assert(n>0);
            
            if this.dim == 3
                s = NaN(this.dim,n);
                phi = 2*pi*rand(1,n);
                s(3,:) = rand(n,1)*2-1;
                r = sqrt(1-s(3,:).^2);
                s(1,:) = r.*cos(phi);
                s(2,:) = r.*sin(phi);
            else
                s=cell2mat(cellfun(@(x){x/norm(x)},num2cell(randn(this.dim,n),1)));
            end
            % Mirror ones with negative on the last dimension up for hemisphere. This
            % may give a disadvantage to ones with exactly zero at the first dimension but
            % since this is hit with quasi probability zero, we neglect
            % this.
            s = (1-2*(s(end,:)<0)).*s;
        end
    end 
end