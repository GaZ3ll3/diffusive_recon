function [curve, sigma_ret, sigma_true] = opt4(sigma_a0)
% k is number of sources
% d is degree of element


global fem Q sigma_s sigma_a u u1 load counter sigma_rate M alpha
% initialize fem 
fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/20000, []', 2);

alpha = 1e-5;

sigma_rate = [];
counter = 1;

% size of problem
n = size(fem.Promoted.nodes, 2);

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

% positive scattering coefficient on element, piecewise constant
sigma_s = 20.0 ;

% initial sigma_a
% sigma_a = 0.1 + 0.05 * (nodes(1,:) + nodes(2, :))';
sigma_a = (0.1 * (1.0 + 0.05 .* (nodes(1, :) > 0.5) .* (nodes(2, :) > 0.5) .* (nodes(1,:) < 0.75) .* (nodes(2, :) < 0.75) +...
    0.15 .* (nodes(1,:) < 0.4) .*(nodes(2, :) < 0.4) .* (nodes(1,:) > 0.2) .* (nodes(2, :) > 0.2)))';
M = fem.assema(1);

if nargin < 1
    sigma_a0 = 0. .* (1 + 2 * (rand(size(sigma_a)) - 0.5));
end
% sigma_a0 = 0.06* ones(n, 1);

% initialize load vector, multiple sources.
source_fcn = @ring;
source_q_nodes = source_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
load = fem.asseml(source_q_nodes);


sigma_t = sigma_a + sigma_s;


qsigma_t = focus_mapping(sigma_t, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a = focus_mapping(sigma_a, fem.Promoted.elems, fem.Facet.Ref');

DSA = fem.assems((1/3)./qsigma_t) + fem.assema(qsigma_a) + 0.5 * Q;
u = DSA\load;

u1 = u .* sigma_a;

% add some noise
u1 = u1 .* (1 + 0.05 * 2 * (rand(size(u)) - 0.5));

% options = optimoptions('fminunc','Display','off','Algorithm',...
%     'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 3000, 'TolFun',...
%     1e-12, 'TolX',1e-16,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');

% problem.options = options;
% problem.x0 = sigma_a0;
% problem.objective = @data;
% problem.solver = 'fminunc';
% 
% 
% 
% sigma_ret = fminunc(problem);

lb  = zeros(size(sigma_a0));  % Lower bound on the variables.
ub  = ones(size(sigma_a0));  % Upper bound on the variables.


sigma_ret = lbfgsb(sigma_a0,lb,ub,'data4','data_grad4',...
           [],'callback','maxiter',1e4,'m',8,'factr',1e-12,...
           'pgtol',1e-12);

fprintf('Relative L2 error of sigma_a is %6.2f \n', norm(sigma_ret - sigma_a)/norm(sigma_a));

sigma_t = sigma_ret + sigma_s;

qsigma_t = focus_mapping(sigma_t, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a = focus_mapping(sigma_ret, fem.Promoted.elems, fem.Facet.Ref');

T = fem.assems((1/3)./qsigma_t) + fem.assema(qsigma_a) + 0.5 * Q;
v = T\load;

fprintf('Relateive L2 error of u is %6.2f \n', norm(u - v)/norm(u));

figure(1);
plot(sigma_rate);
curve = sigma_rate;
sigma_true = sigma_a;

figure(2);
trimeshC1(sigma_a);

figure(3);
trimeshC1(sigma_ret);

figure(4);
trimeshC1(sigma_ret - sigma_a);

end



function trimeshC1(sigma)

global fem

numofnodes = size(fem.Promoted.nodes, 2);
trimesh(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), sigma(1:numofnodes));colorbar;colormap jet;

end
