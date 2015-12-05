function [curve, sigma_ret, sigma_true] = opt5(sigma_a0)
% k is number of sources
% d is degree of element


global fem Q sigma_s sigma_a Gamma u_1 u_2 u1 u2 load_1 load_2 counter sigma_rate gamma_rate M alpha beta_
% initialize fem 
fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/4000, []', 2);

alpha = 1e-6;
beta_ = 1e-6;

sigma_rate = [];
gamma_rate = [];
counter = 1;

% size of problem
% n = size(fem.Promoted.nodes, 2);

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
% sigma_a = 0.5 + 0.2 * (nodes(1,:) - nodes(2, :))';
sigma_a = (0.1 * (1.0 + 0.05 .* (nodes(1, :) > 0.5) .* (nodes(2, :) > 0.5) .* (nodes(1,:) < 0.75) .* (nodes(2, :) < 0.75) +...
    0.15 .* (nodes(1,:) < 0.4) .*(nodes(2, :) < 0.4) .* (nodes(1,:) > 0.2) .* (nodes(2, :) > 0.2)))';
M = fem.assema(1);

% Gamma = 1.0  + 0.2 * (nodes(2, :) - nodes(1,:))';
Gamma = 0.5 * ones(size(sigma_a, 1), 1) +...
    (0.05 .* (nodes(1, :) > 0.5) .* (nodes(2, :) > 0.5) .* (nodes(1,:) < 0.75) .* (nodes(2, :) < 0.75))';

if nargin < 1
    sigma_a0 = sigma_a.* (1 + 0.2 * (rand(size(sigma_a)) - 0.5));
    Gamma_0 = Gamma .* (1 + 0.2 * (rand(size(sigma_a)) - 0.5));
end
% sigma_a0 = 0.06* ones(n, 1);

% initialize load vector, multiple sources.
source_fcn_1 = @SinSource;
source_q_nodes_1 = source_fcn_1(fem.Qnodes(1,:), fem.Qnodes(2,:));
load_1 = fem.asseml(source_q_nodes_1);

source_fcn_2 = @CosSource;
source_q_nodes_2 = source_fcn_2(fem.Qnodes(1,:), fem.Qnodes(2,:));
load_2 = fem.asseml(source_q_nodes_2);


sigma_t = sigma_a + sigma_s;


qsigma_t = focus_mapping(sigma_t, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a = focus_mapping(sigma_a, fem.Promoted.elems, fem.Facet.Ref');

DSA = fem.assems((1/3)./qsigma_t) + fem.assema(qsigma_a) + 0.5 * Q;
u_1 = DSA\load_1;
u_2 = DSA\load_2;

u1 = u_1 .* sigma_a .* Gamma;
u2 = u_2 .* sigma_a .* Gamma;
% add some noise
% u1 = u1 .* (1 + 0.05 * 2 * (rand(size(u1)) - 0.5));
% u2 = u2 .* (1 + 0.05 * 2 * (rand(size(u2)) - 0.5));
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

lb  = zeros(size(sigma_a0, 1) * 2, 1);  % Lower bound on the variables.
ub  = 10 * ones(size(sigma_a0, 1) * 2, 1);  % Upper bound on the variables.

start = [sigma_a0; Gamma_0];


ret = lbfgsb(start,lb,ub,'data5','data_grad5',...
           [],'callback5','maxiter',1e4,'m',8,'factr',1e-12,...
           'pgtol',1e-12);

       
sigma_ret = ret(1:size(u1, 1));
Gamma_ret = ret(size(u1, 1) + 1:end);
       
fprintf('Relative L2 error are %6.2f, %6.2f \n',...
    norm(sigma_ret - sigma_a)/norm(sigma_a), ...
    norm(Gamma_ret - Gamma)/norm(Gamma));

sigma_t = sigma_ret + sigma_s;

qsigma_t = focus_mapping(sigma_t, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a = focus_mapping(sigma_ret, fem.Promoted.elems, fem.Facet.Ref');

T = fem.assems((1/3)./qsigma_t) + fem.assema(qsigma_a) + 0.5 * Q;
v1 = T\load_1;
v2 = T\load_2;

v1 = sigma_ret .* v1 .* Gamma_ret;
v2 = sigma_ret .* v2 .* Gamma_ret;

fprintf('Relateive L2 error of u are %6.2f, %6.2f \n',...
    norm(u1 - v1, inf)/norm(u1, inf),...
    norm(u2 - v2, inf)/norm(u2, inf));

figure(1);
plot(sigma_rate);
figure(8);
plot(gamma_rate);
curve = sigma_rate;
sigma_true = sigma_a;

figure(2);
trimeshC1(sigma_a);

figure(3);
trimeshC1(sigma_ret);

figure(4);
trimeshC1(sigma_ret - sigma_a);

figure(5);
trimeshC1(Gamma);

figure(6);
trimeshC1(Gamma_ret);

figure(7);
trimeshC1(Gamma - Gamma_ret);

figure(9);
trimeshC1(Gamma.*sigma_a - Gamma_ret .* sigma_ret);

end



function trimeshC1(sigma)

global fem

numofnodes = size(fem.Promoted.nodes, 2);
trimesh(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), sigma(1:numofnodes));colorbar;colormap jet;

end

function val = SinSource(x, y)
    val = 10 * (0.2 + sin(pi/4 * ((x - 0.5).^2 + (y - 0.5).^2)));
end

function val = CosSource(x, y)
    val = 10 * cos(pi/4 * ((x - 0.5).^2 + (y - 0.5).^2));
end
