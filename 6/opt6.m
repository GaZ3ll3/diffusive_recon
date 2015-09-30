classdef opt6 < handle
    
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
    end
    
    properties(Access = private)
        sigma_a_local
        gamma_local
        
        data_local
        sols_local
        

    end
    
    
    methods
        function this = opt6()
            % constructor
            
            % finite element with out pml.
            this.fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/2000, []', 2);
            
            boundary = Boundary();
            boundary.set_boundary('x - 1');
            boundary.set_boundary('y - 1');
            boundary.set_boundary('x');
            boundary.set_boundary('y');

            [bc1, bc2, bc3, bc4] = boundary.get_boundary(this.fem.Promoted.edges, this.fem.Promoted.nodes, 4);

            boundary.setDirichlet(bc1);
            boundary.setDirichlet(bc2);
            boundary.setDirichlet(bc3);
            boundary.setDirichlet(bc4);

            % boundary integral term, constant
            this.Edge = this.fem.assemlbc(1, bc1) + this.fem.assemlbc(1, bc2) + this.fem.assemlbc(1, bc3) + this.fem.assemlbc(1, bc4);
            % matrix
            this.Mass = this.fem.assema(1);
            this.Stiff = this.fem.assems(1);
            
%             this.source{1} = @this.SinSource;
%             this.source{2} = @this.CosSource;
            
            this.alpha{1} = 1e-5;
            this.alpha{2} = 1e-5;
            
            this.rate = {[], []};
            
            nodes = this.fem.Promoted.nodes;
            % variable
            this.sigma_a = this.SigmaFcn(nodes(1,:), nodes(2,:));
            this.gamma   = this.GammaFcn(nodes(1,:), nodes(2,:));
            this.sigma_s = 5.0;

            % initials
            this.sigma_a_0 = this.sigma_a .* ( 1 + 0.5 * (rand(size(this.sigma_a)) - 0.5));
            this.gamma_0   = this.gamma .* (1 + 0.5 * (rand(size(this.gamma)) - 0.5));
%             this.sigma_a_0 = 0.1 * ones(size(this.sigma_a));
%             this.gamma_0   = 0.5 * ones(size(this.gamma));
            
%             this.load{1} = this.fem.asseml(this.source{1}(this.fem.Qnodes(1, :), this.fem.Qnodes(2, :)));
%             this.load{2} = this.fem.asseml(this.source{2}(this.fem.Qnodes(2, :), this.fem.Qnodes(2, :)));
            for i = 1: 2
                this.load{i} = rand(size(this.sigma_a_0));
            end
            
            sigma_t = this.sigma_a + this.sigma_s;
            qsigma_t = this.mapping(sigma_t, this.fem.Promoted.elems, this.fem.Facet.Ref');
            qsigma_a = this.mapping(this.sigma_a, this.fem.Promoted.elems, this.fem.Facet.Ref');
            
            DSA = this.fem.assems((1/3)./qsigma_t) + this.fem.assema(qsigma_a) + 0.5 * this.Edge;
            
%             this.sols{1} = DSA\this.load{1};
%             this.sols{2} = DSA\this.load{2};

            for i = 1 :2
                this.sols{i} = DSA\this.load{i};
                this.data{i} = this.sols{i} .* this.gamma .* this.sigma_a;
            end
            
            
%             this.data{1} = this.sols{1} .* this.gamma .* this.sigma_a;
%             this.data{2} = this.sols{2} .* this.gamma .* this.sigma_a;
            
        end
        
        function delete(this)
            % override destructor
            this.fem.delete();
        end
        
        function [f, g] = objective_gradient(this, ret)
            this.sigma_a_local = ret(1:end/2);
            this.gamma_local   = ret(end/2 + 1: end);
            
            sigma_t = this.sigma_a_local + this.sigma_s;
            
            qsigma_t = this.mapping(sigma_t, this.fem.Promoted.elems, this.fem.Facet.Ref');
            qsigma_a = this.mapping(this.sigma_a_local, this.fem.Promoted.elems, this.fem.Facet.Ref');
            
            DSA = this.fem.assems((1/3)./qsigma_t) + this.fem.assema(qsigma_a) + 0.5 * this.Edge;
            
            tmp = cell(2, 1);
            
            for i = 1 : 2
                this.sols_local{i} = DSA\this.load{i};
                this.data_local{i} = this.sols_local{i} .* this.sigma_a_local .* this.gamma_local;
                tmp{i} = DSA\(this.sigma_a_local .* (this.Mass * (this.data{i} - this.data_local{i})));
            end
                        
            r{1} = 0.5 * this.sigma_a_local' * this.Stiff * this.sigma_a_local;
            h{1} = this.Stiff * this.sigma_a_local;
            r{2} = 0.5 * this.gamma_local' * this.Stiff * this.gamma_local;
            h{2} = this.Stiff * this.gamma_local;
            

            

            
            g{1} = zeros(size(this.sigma_a));
            g{2} = zeros(size(this.gamma));
            
            for i = 1 : 2
                g{1} = g{1} + ...
                    this.fem.assemnode(tmp{i}, this.sols{i}, -1/3./sigma_t./sigma_t, ones(size(sigma_t))) - ... 
                this.sols_local{i} .* (this.Mass * (this.data{i} - this.data_local{i}));
            
                g{2} = g{2} - (this.sigma_a_local .* this.sols{i}) .* (this.Mass * (this.data{i} - this.data_local{i}));
            end
            
            g{1} = this.gamma_local .* ...
                (g{1}) + this.alpha{1} * h{1}; 
            g{2} = g{2} +...
                this.alpha{2} * h{2};
            
            g = [g{1}; g{2}];  

            f = 0.5 * (this.data{1} - this.data_local{1})' * this.Mass * (this.data{1} - this.data_local{1}) + ...
                0.5 * (this.data{2} - this.data_local{2})' * this.Mass * (this.data{2} - this.data_local{2}) + ...
                this.alpha{1} * r{1} + this.alpha{2} * r{2};
            
            this.rate{1} = [this.rate{1}, norm(this.sigma_a - this.sigma_a_local)/norm(this.sigma_a)];
            this.rate{2} = [this.rate{2}, norm(this.gamma - this.gamma_local)/norm(this.gamma)];
        end
        

        
        function optimize(this,start)
            
            if (nargin == 1)
                start = [this.sigma_a_0; this.gamma_0];
            end
             
            options = optimoptions('fminunc','Display','iter','Algorithm',...
            'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 3000, 'TolFun',...
            1e-12, 'TolX',1e-16,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');

            problem.options = options;
            problem.x0 = start;
            problem.objective = @this.objective_gradient;
            problem.solver = 'fminunc';

            ret = fminunc(problem);
            
            this.sigma_a_ = ret(1:end/2);
            this.gamma_   = ret(end/2 + 1: end);
        end
        
        function plot(this, level)
            
            if (level > 1)
            
            figure(1);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
        
            end
            
            if (level > 2)
            figure(3);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(4);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            end
            
            if (level > 1)
            figure(5);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes) - this.sigma_a_(1:numofnodes)...
            ,'EdgeColor', 'None');shading interp;colorbar;colormap jet;
            figure(6);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.gamma(1:numofnodes) -  this.gamma_(1:numofnodes)...
            ,'EdgeColor', 'None');shading interp;colorbar;colormap jet;
        
            end
            
            if (level > 0)
            figure(7);
            plot(this.rate{1});
            figure(8);
            plot(this.rate{2});
            end
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
            val = 0.1 +0.05 * (sin(4 * pi *  x) .* sin(4 * pi * y))';
        end
        %gamma
        function val = GammaFcn(x, y)
            val = 0.5  + 0.2 * (sin(4 * pi * (x + 0.25)) .* sin(4 * pi * (y)))';         
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
    end
    
end

