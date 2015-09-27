function [f] = data5(ret)

global fem u1 u2 load_1 load_2 sigma_s sigma_a Gamma Q counter sigma_rate gamma_rate M alpha beta_


sigma_a_ = ret(1: size(u1, 1));
Gamma_ret = ret(size(u1, 1) + 1: end);

sigma_t_ = sigma_a_ + sigma_s;

qsigma_t_ = focus_mapping(sigma_t_, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a_ = focus_mapping(sigma_a_, fem.Promoted.elems, fem.Facet.Ref');

T = fem.assems((1/3)./qsigma_t_) + fem.assema(qsigma_a_) + 0.5 * Q;

v1 = T\load_1;
v2 = T\load_2;

v1 = v1 .* sigma_a_ .* Gamma_ret;
v2 = v2 .* sigma_a_ .* Gamma_ret;

[r1, ~] = reg5(sigma_a_);
[r2, ~] = reg5(Gamma_ret);


f = 0.5 * (u1 - v1)' * M *  (u1 - v1) +...
    0.5 *(u2 - v2)' * M * (u2 - v2) +....
    alpha * (r1) + beta_ * r2;


sigma_rate(counter) = norm(sigma_a_ - sigma_a)/norm(sigma_a);
gamma_rate(counter) = norm(Gamma_ret - Gamma)/norm(Gamma);
counter = counter + 1;


end