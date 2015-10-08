classdef radiative2 < handle
    %OPT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        fem 
        dom
        
        % now these are structs
        sol
        sol_
        source
        
        data
        data_
        
        Mass
        Stiff
        alpha
        
        m_
        sigma_t_
        Q_
        Qt_
        
        
        
        sigma_a
        sigma_s
    
        sigma_a_0
        sigma_s_0
        sigma_s_
        sigma_a_
    end
    
    
    
    methods 
        function this = radiative2()
            this.fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 64 * 64), []');
            this.dom = DOM(16);
            
            this.dom.rayint(this.fem.Promoted.nodes, this.fem.Promoted.elems, this.fem.Promoted.neighbors);
            % two instances of sources
            this.source{1} = this.source_fcn1(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            this.source{2} = this.source_fcn2(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            
            
            
            this.sigma_a = this.Sigma_a_Fcn(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            this.sigma_s = this.Sigma_s_Fcn(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            
            this.alpha = 1e-5;
            this.Mass = this.fem.assema(1);
            this.Stiff =  this.fem.assems(1);
            
%             this.sigma_a_0 = this.sigma_a .* (1.0 + 0.8* (rand(size(this.sigma_s)) - 0.5));
%             this.sigma_s_0 = 5.0 * ones(size(this.sigma_s));
            this.sigma_a_0 = 0.1 * ones(size(this.sigma_a));
            
            sigma_t = this.sigma_a + this.sigma_s;
            
            m = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, sigma_t);
            
            
            % two instances, sensitive ?
            [this.sol{1},~,~,~,~]= gmres(speye(size(m, 1)) - m' * sparse(1:size(m,1), 1:size(m,1), this.sigma_s), m' * (this.source{1}),10 , 1e-14, 400);
            [this.sol{2},~,~,~,~]= gmres(speye(size(m, 1)) - m' * sparse(1:size(m,1), 1:size(m,1), this.sigma_s), m' * (this.source{2}),10 , 1e-14, 400);
            this.data = this.sol{1}./this.sol{2};
            
        end
        
        function f = objective(this, ret) 
            
            this.sigma_t_ = this.sigma_s + ret;
            this.m_ = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, this.sigma_t_);
            this.Q_ = speye(size(this.m_, 1)) - this.m_' * sparse(1:size(this.m_, 1), 1:size(this.m_, 1), this.sigma_s);
            this.Qt_ = speye(size(this.m_, 1)) - this.m_ * sparse(1:size(this.m_, 1), 1:size(this.m_, 1), this.sigma_s);
            
            [this.sol_{1}, ~, ~, ~, ~] = gmres(this.Q_, this.m_' * (this.source{1}), 10, 1e-14, 400);
            [this.sol_{2}, ~, ~, ~, ~] = gmres(this.Q_, this.m_' * (this.source{2}), 10, 1e-14, 400);
            
            this.data_ = this.sol_{1} ./this.sol_{2};

               
            r = 0.5 * ret' * this.Stiff * ret;
            f = 0.5 * (this.data_ - this.data)' * this.Mass * (this.data_ - this.data) + this.alpha * r;
            
        end
        
        function g = gradient(this, ret)
            
            
            
            
            
            
            tmp = this.Mass * (this.data_ - this.data);
            
            tmp1 = tmp./this.sol_{2};
            tmp2 = tmp.*this.data_./this.sol_{2};
            
            h = this.Stiff * ret;
            [L1, ~, ~,~,~] = gmres(this.Qt_, tmp1, 10, 1e-14, 400); 
            [L2, ~, ~,~,~] = gmres(this.Qt_, tmp2, 10, 1e-14, 400); 
            
            g1 = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
                ...
                 L1,(this.source{1} + this.sigma_s .* this.sol_{1}), this.sigma_t_);
             
            g2 = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
                ...
                 L2,(this.source{2} + this.sigma_s .* this.sol_{2}), this.sigma_t_);
             
            g = g1 - g2 + this.alpha * h;
             
             
             
             
             
             
             
             
             
             
             
        end
        
        function [f, g] = objective_gradient(this, ret)
            
            sigma_t = this.sigma_s + ret;
            m = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, sigma_t);
            
            Q = speye(size(m, 1)) - m' * sparse(1:size(m, 1), 1:size(m, 1), this.sigma_s);
            Qt = speye(size(m, 1)) - m * sparse(1:size(m, 1), 1:size(m, 1), this.sigma_s);
            
            
            this.sol_ = (Q)\(m' * (this.source));
            

            r = 0.5 * ret' * this.Stiff * ret;
            h = this.Stiff * ret;
            
            f = 0.5 * (this.sol_ - this.sol)' * this.Mass * (this.sol_ - this.sol) + this.alpha * r;
            tmp = this.Mass * (this.sol_ - this.sol);
             
            L = Qt\tmp; 
            
            g = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
                ...
                L,(this.source + this.sigma_s .* this.sol_),  sigma_t) + this.alpha *h;
        end
        
        

        function LBFGS(this, objstr, gradstr, cbstr, start)
            if (nargin == 4)
                start = this.sigma_a_0;
            end
            lb = 0. * ones(size(start));
            ub = Inf * ones(size(start));
            
            ret = lbfgsb(start,lb,ub, objstr, gradstr,...
                    [],cbstr,'maxiter',1e4,'m',50,'factr',1e-16,...
                    'pgtol',1e-16);
       
            this.sigma_a_ = ret;            
        end    

        
        function optimize(this, start)
            
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
            
            
        end
        
        function plot(this)
            figure(1);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(3);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(4);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
        end
        
        function output(this)   
            h = figure(1);set(h, 'Position', [200, 200, 1350, 700]);
            subplot(2, 3, 1);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;view(2);title('reconstructed $\tilde{\sigma}_a$','Interpreter','latex');
            subplot(2, 3, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;view(2);title('true $\sigma_a$','Interpreter','latex');
            subplot(2, 3, 3);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes) - this.sigma_a_(1:numofnodes),...
            'EdgeColor', 'None');shading interp;title('error of $\sigma_a - \tilde{\sigma}_a$','Interpreter','latex')
            colorbar;colormap jet;view(2);
            subplot(2, 3, 4);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;view(2);title('solution from $\tilde{\sigma}_a$', 'Interpreter', 'latex');
            subplot(2, 3, 5);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;view(2);title('solution from ${\sigma}_a$', 'Interpreter', 'latex');
            subplot(2, 3, 6);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol_(1:numofnodes) - this.sol(1:numofnodes),...
            'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;view(2);title('error of solutions', 'Interpreter', 'latex');
            
        end
        
    end 
    
    methods (Static)
        
        function val = source_fcn1(x, y)
            val = 10 + x + y;
        end
        
        function val = source_fcn2(x, y)
           val = 5 + 8 * ((x - 0.5) .^ 2 + (y - 0.5) .^2); 
        end
        
        function val = Sigma_a_Fcn(x, y)
%             val =(0.1 * (1.0 + 0.45 .* (x > 0.5)...
%                 .* (y > 0.5) .* (x < 0.75)...
%                 .* (y < 0.75) +...
%                 0.45 .* (x < 0.4) .*(y < 0.4)...
%                 .* (x > 0.2) .* (y > 0.2)));  
            val = 0.1 * (2 + sin(4 * pi * x) .* sin(4 * pi * y));
%             val = 0.1 * ones(size(x));
        end
        
        function val = Sigma_s_Fcn(x, y)
            val = 5.0 * ones(size(x));
        end
        
    end
end