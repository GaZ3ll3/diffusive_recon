function [curve, sigma_ret, sigma_true] = opt(sigma_a0)
% k is number of sources
% d is degree of element


global fem Q sigma_s sigma_a u load counter sigma_rate
% initialize fem 
fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/5000, []', 2);

sigma_rate = [];
counter = 1;

% size of problem
n = size(fem.Promoted.elems, 2);
% m = size(fem.Promoted.nodes, 2);

% initiaize boundary condition
boundary = Boundary();
boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');

[bc1, bc2, bc3, bc4] = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);

boundary.setDirichlet(bc1);
boundary.setDirichlet(bc2);
boundary.setDirichlet(bc3);
boundary.setDirichlet(bc4);

% boundary integral term, constant
Q = fem.assemlbc(1, bc1) + fem.assemlbc(1, bc2) + fem.assemlbc(1, bc3) + fem.assemlbc(1, bc4);

nodes = fem.Promoted.nodes;

X = (nodes(1, fem.Promoted.elems(1, :))' + nodes(1, fem.Promoted.elems(2, :))' +nodes(1, fem.Promoted.elems(3, :))')/3.0 ;
Y = (nodes(2, fem.Promoted.elems(1, :))' + nodes(2, fem.Promoted.elems(2, :))' + nodes(2, fem.Promoted.elems(3, :))')/3.0;

% positive scattering coefficient on element, piecewise constant
sigma_s = 20.0 ;

% initial sigma_a
sigma_a = 0.1 * (1.0 + 0.05 .* (X > 0.5) .* (Y > 0.5) .* (X < 0.75) .* (Y < 0.75) +...
    0.15 .* (X < 0.4) .*(Y < 0.4) .* (X > 0.2) .* (Y > 0.2));


if nargin < 1
    sigma_a0 = 0.1 * zeros(n, 1);
end
% sigma_a0 = 0.06* ones(n, 1);

% initialize load vector, multiple sources.
source_fcn = @ring;
source_q_nodes = source_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
load = fem.asseml(source_q_nodes);


sigma_t = sigma_a + sigma_s;
DSA = fem.assems((1/3)./sigma_t) + fem.assema(sigma_a) + 0.5 * Q;
u = DSA\load;


% add some noise
u = u .* (1 + 0.01 * 2 * (rand(size(u)) - 0.5));

options = optimoptions('fminunc','Display','iter','Algorithm',...
    'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 100, 'TolFun',...
    1e-12, 'TolX',1e-16,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');

problem.options = options;
problem.x0 = sigma_a0;
problem.objective = @data;
problem.solver = 'fminunc';



sigma_ret = fminunc(problem);

fprintf('Relative L2 error of sigma_a is %6.2f percentage \n', norm(sigma_ret - sigma_a)/norm(sigma_a) * 100);

sigma_t = sigma_ret + sigma_s;
T = fem.assems((1/3)./sigma_t) + fem.assema(sigma_a) + 0.5 * Q;
v = T\load;

fprintf('Relateive L2 error of u is %6.2f percentage \n', norm(u - v)/norm(u));

figure(1);
plot(sigma_rate);
curve = sigma_rate;
sigma_true = sigma_a;

figure(2);
trimeshC0(sigma_a);

figure(3);
trimeshC0(sigma_ret);

end


function [f, g] = data(sigma_a_)

global fem u load sigma_s sigma_a Q counter sigma_rate

sigma_t_ = sigma_a_ + sigma_s;

T = fem.assems((1/3)./sigma_t_) + fem.assema(sigma_a_) + 0.5 * Q;

M = fem.assema(1);

v = T\load;

alpha = 1e-4;

[r, h] = reg(sigma_a_);

f = 0.5 * (u - v)' * M *  (u - v) + alpha * r;

sigma_rate(counter) = norm(sigma_a_ - sigma_a)/norm(sigma_a);
counter = counter + 1;

tmp = T\(M * (u - v));

g = fem.assemelem(tmp, v, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + alpha * h;


end


function [r, g] = reg(sigma_a_)

global fem


r = 0.;
n = size(fem.Promoted.neighbors, 2);
g = zeros(n , 1);
for i = 1:n
    for j = 1:3
        t = fem.Promoted.neighbors(j, i);
        if (t ~= 0 && t > i)
            r = r + (sigma_a_(i) - sigma_a_(t))^2;
            g(i) = g(i) + 2 * (sigma_a_(i) -sigma_a_(t));
            g(t) = g(t) - 2 * (sigma_a_(i) - sigma_a_(t));
        end
    end
end


end

function trimeshC0(sigma)

global fem

F = fem.TriMesh';
V = fem.Promoted.nodes(1:2, :)';
C = V(F(:),:);
SC = repmat(sigma, size(F, 2), 1);
FC = bsxfun(@plus,size(F,1)*(0:size(F,2)-1),(1:size(F,1))');
trisurf(FC, C(:,1), C(:,2), SC,'EdgeColor','None');colorbar;colormap jet;

end
