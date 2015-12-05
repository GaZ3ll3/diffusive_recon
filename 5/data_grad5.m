function [g] =  data_grad5(ret)
global fem u1 u2 u_1 u_2 load_1 load_2 sigma_s Q M alpha beta_


sigma_a_ = ret(1:size(u1, 1));
Gamma_ret = ret(size(u1, 1) + 1: end);

sigma_t_ = sigma_a_ + sigma_s;

qsigma_t_ = focus_mapping(sigma_t_, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a_ = focus_mapping(sigma_a_, fem.Promoted.elems, fem.Facet.Ref');

T = fem.assems((1/3)./qsigma_t_) + fem.assema(qsigma_a_) + 0.5 * Q;
v1 = T\load_1;
v2 = T\load_2;

mv1 = sigma_a_ .* v1 .* Gamma_ret;
mv2 = sigma_a_ .* v2 .* Gamma_ret;

tmp1 = T\(sigma_a_ .* (M * (u1 - mv1)));
tmp2 = T\(sigma_a_ .* (M * (u2 - mv2)));

[~, h1] = reg5(sigma_a_);
[~, h2] = reg5(Gamma_ret);

g1 = Gamma_ret .* ...
    (fem.assemnode(tmp1, v1, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + ...
    -v1 .* (M * (u1 - mv1)) + ...
    fem.assemnode(tmp2, v2, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + ...
    -v2 .* (M * (u2 - mv2)) )+ alpha * h1;

g2 = -(sigma_a_ .* u_1) .* (M *(u1 - mv1)) -(sigma_a_ .* u_2) .* (M * (u2 - mv2)) + beta_* h2;

g = [g1; g2];

end