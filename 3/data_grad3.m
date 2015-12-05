function [g] =  data_grad3(sigma_a_)
global fem u load sigma_s Q M alpha

sigma_t_ = sigma_a_ + sigma_s;

qsigma_t_ = focus_mapping(sigma_t_, fem.Promoted.elems, fem.Facet.Ref');
qsigma_a_ = focus_mapping(sigma_a_, fem.Promoted.elems, fem.Facet.Ref');

T = fem.assems((1/3)./qsigma_t_) + fem.assema(qsigma_a_) + 0.5 * Q;
v = T\load;

tmp = T\((M * (u - v)));
% tmp = T\(log(u) - log(v));

[~, h] = reg3(sigma_a_);

g = fem.assemnode(tmp, v, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + alpha * h;
% g = fem.assemnode(tmp, v, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + alpha * h;

end