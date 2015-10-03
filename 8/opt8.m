classdef opt8 < handle
    %OPT8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        fem
        
        sols
        data
        source
        load
        
        % variable
        sigma_a
        sigma_s
        gamma
        
        % initial variable
        sigma_a_0
        gamma_0
        
        
        % output vaiable
        sigma_a_
        gamma_
        
        %matrix
        Edge
        Mass
        Stiff
        
        % optimization
        rate
        alpha
        
        
        % quotient
        quotient
    end
    
    properties(Access = private)
        sigma_a_local
        gamma_local
        
        data_local
        sols_local
        
        quotient_local
    end
    
    methods
        function this = opt8()
            % constructor
            
            % finite element with out pml.
            this.fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/4000, []', 3);
            
            boundary = Boundary();
            boundary.set_boundary('x - 1');
            boundary.set_boundary('y - 1');
            boundary.set_boundary('x');
            boundary.set_boundary('y');

            [bc1, bc2, bc3, bc4] = boundary.get_boundary(...
                this.fem.Promoted.edges, this.fem.Promoted.nodes, 4);

            boundary.setDirichlet(bc1);
            boundary.setDirichlet(bc2);
            boundary.setDirichlet(bc3);
            boundary.setDirichlet(bc4);

            % boundary integral term, constant
            this.Edge = this.fem.assemlbc(1, bc1) + this.fem.assemlbc(1, bc2) +...
                this.fem.assemlbc(1, bc3) + this.fem.assemlbc(1, bc4);
            % matrix
            this.Mass = this.fem.assema(1);
            this.Stiff = this.fem.assems(1);
            
            this.source{1} = @this.SinSource;
            this.source{2} = @this.CosSource;
            
            this.alpha{1} = 5e-9;
            this.alpha{2} = 5e-9;
            
            this.rate = {[], []};
            
            nodes = this.fem.Promoted.nodes;
            % variable
            this.sigma_a = this.SigmaFcn(nodes(1,:), nodes(2,:));
            this.gamma   = this.GammaFcn(nodes(1,:), nodes(2,:));
            this.sigma_s = 20.0;

            % initials
            %this.sigma_a_0 = this.sigma_a .* ( 1 + 0.4 * (rand(size(this.sigma_a)) - 0.5));
            %this.sigma_a_0 = 0.1 * ones(size(this.sigma_a));
            this.sigma_a_0 = this.sigma_a;
            this.gamma_0   = this.gamma .* (1 + 0.2 * (rand(size(this.gamma)) - 0.5));
            
            
            this.load{1} = this.fem.asseml(this.source{1}(this.fem.Qnodes(1, :), this.fem.Qnodes(2, :)));
            this.load{2} = this.fem.asseml(this.source{2}(this.fem.Qnodes(2, :), this.fem.Qnodes(2, :)));
            
            sigma_t = this.sigma_a + this.sigma_s;
            qsigma_t = this.mapping(sigma_t, this.fem.Promoted.elems, this.fem.Facet.Ref');
            qsigma_a = this.mapping(this.sigma_a, this.fem.Promoted.elems, this.fem.Facet.Ref');
            
            DSA = this.fem.assems((1/3)./qsigma_t) + this.fem.assema(qsigma_a) + 0.5 * this.Edge;
            
            this.sols{1} = DSA\this.load{1};
            this.sols{2} = DSA\this.load{2};
            
            this.data{1} = this.sols{1} .* this.gamma .* this.sigma_a;
            this.data{2} = this.sols{2} .* this.gamma .* this.sigma_a;
            
            
            % data for reconstructing sigma_a
            this.quotient = this.data{1}.*(1 + 0.001 * (rand(size(this.sigma_a))- 0.5)) ./...
                this.data{2};
            
        end
        
        function [f, g] = objective_gradient(this, ret)
            % extract information
            this.sigma_a_local = ret;
            %this.gamma_local   = ret(end/2 + 1: end);
            
            sigma_t = this.sigma_a_local + this.sigma_s;
            
            qsigma_t = this.mapping(sigma_t, this.fem.Promoted.elems, this.fem.Facet.Ref');
            qsigma_a = this.mapping(this.sigma_a_local, this.fem.Promoted.elems, this.fem.Facet.Ref');
            
            DSA = this.fem.assems((1/3)./qsigma_t) + this.fem.assema(qsigma_a) + 0.5 * this.Edge;
            
            this.sols_local{1} = DSA\this.load{1};
            this.sols_local{2} = DSA\this.load{2};
            
            % quotient data
            this.quotient_local = this.sols_local{1} ./ this.sols_local{2};
            
            carry = this.Mass * (this.quotient - this.quotient_local);
            
            tmp{1} = DSA\( carry .* this.sols_local{1} ./ this.sols_local{2} ./this.sols_local{2});
            tmp{2} = DSA\( carry ./ this.sols_local{2});
                        
            r{1} = 0.5 * this.sigma_a_local' * this.Stiff * this.sigma_a_local;
            h{1} = this.Stiff * this.sigma_a_local;
%             r{2} = 0.5 * this.gamma_local' * this.Stiff * this.gamma_local;
%             h{2} = this.Stiff * this.gamma_local;
            
            g = - this.fem.assemnode(tmp{1}, this.sols{2}, -1/3./sigma_t./sigma_t, ones(size(sigma_t))) +...
                this.fem.assemnode(tmp{2}, this.sols{1}, -1/3./sigma_t./sigma_t, ones(size(sigma_t))) + ...
                this.alpha{1} * h{1};

            f = 0.5 * (this.quotient - this.quotient_local)' * this.Mass * (this.quotient - this.quotient_local) +...
                this.alpha{1} * r{1};
            
            this.rate{1} = [this.rate{1}, norm((this.sigma_a - this.sigma_a_local)./this.sigma_a, inf)];
            this.rate{2} = [this.rate{2}, norm(this.sols_local{1} .* this.sigma_a_local ./(this.sols{1} .* this.sigma_a), inf) - 1];
        end
 
        function optimize(this,start)
            
            if (nargin == 1)
                start = this.sigma_a_0;
            end
            options = optimoptions('fminunc','Display','iter','Algorithm',...
            'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 800, 'TolFun',...
            1e-20, 'TolX',1e-20,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');

            problem.options = options;
            problem.x0 = start;
            problem.objective = @this.objective_gradient;
            problem.solver = 'fminunc';

            ret = fminunc(problem);
            
            this.sigma_a_ = ret;
        %    this.gamma_   = ret(end/2 + 1: end);
        end
        
        function LBFGS(this, objstr, gradstr, cbstr)
            start = this.sigma_a_0;
            lb = -5 * ones(size(start));
            ub = 5 * ones(size(start));
            
            ret = lbfgsb(start,lb,ub, objstr, gradstr,...
                    [],cbstr,'maxiter',1e5,'m',12,'factr',1e-12,...
                    'pgtol',1e-16);
       
            this.sigma_a_ = ret;
           % this.gamma_   = ret(end/2 + 1: end);
            
        end        
        
        
        function delete(this)
            % override destructor
            this.fem.delete();
        end
        
        
        function plot(this)
            
            figure(1);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
%             figure(2);
%             trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
%             this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma_(1:numofnodes),'EdgeColor', 'None');shading interp;
%             colorbar;colormap jet;
        
            figure(3);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
%             figure(4);
%             trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
%             this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma(1:numofnodes),'EdgeColor', 'None');shading interp;
%             colorbar;colormap jet;
        
            figure(5);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes) - this.sigma_a_(1:numofnodes)...
            ,'EdgeColor', 'None');shading interp;colorbar;colormap jet;
%             figure(6);
%             trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
%             this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma(1:numofnodes) -  this.gamma_(1:numofnodes)...
%             ,'EdgeColor', 'None');shading interp;colorbar;colormap jet;
        
            figure(7);
            plot(this.rate{1});
            figure(8);
            plot(this.rate{2});
        end
    end
    
    methods(Static)
        % all methods should return column vector
        % source
        function val = SinSource(x, y)
            val = 10 * (0.2 + sin(pi/4 * ((x - 0.5).^2 + (y - 0.5).^2)));
        end
        % source
        function val = CosSource(x, y)
            val = 10 * cos(pi/4 * ((x - 0.5).^2 + (y - 0.5).^2));
        end
        % sigma_a
        function val = SigmaFcn(x, y)
            val = (0.2 * (1.0 + 0.25 .* (x > 0.5)...
                .* (y > 0.5) .* (x < 0.75)...
                .* (y < 0.75) +...
                0.45 .* (x < 0.4) .*(y < 0.4)...
                .* (x > 0.2) .* (y > 0.2)))';
%             val = 0.1 * (1.0 + 0.15 * (sin(4 * pi * x) .* sin(4 * pi * y)))';
        end
        %gamma
        function val = GammaFcn(x, y)
            val = 0.5  +...
                (0.05 .* (x > 0.5) .* ...
                (y > 0.5) .* (x < 0.75) ...
                .* (y < 0.75))';
%             val = 0.5 * (1.0 + 0.2 * (sin(4 * pi * (x + 0.25)) .* sin(4 * pi * y)))';
        end
        % interpolation
        function [interpolate] = mapping(func, elems, trans_ref)

            % allocate memory
            numberofqnodes = size(trans_ref, 1);

            interpolate = zeros(numberofqnodes, size(elems, 2));

            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        
        % callback
        function callback(t, f, x)
            fprintf('%5d\t\t%0.8g\n', t, f);
        end
    end
    
end

