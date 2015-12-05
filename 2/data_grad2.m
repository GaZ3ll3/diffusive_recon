function [g] =  data_grad2(sigma_a_)
global fem u load sigma_s Q M alpha

sigma_t_ = sigma_a_ + sigma_s;

T = fem.assems((1/3)./sigma_t_) + fem.assema(sigma_a_) + 0.5 * Q;
v = T\load;
% tmp = T\(M * (u - v));
tmp = T\(log(u) - log(v));
[~, h] = reg2(sigma_a_);
% g = fem.assemelem(tmp, v, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + alpha * h;
g = fem.assemelem(tmp, v, -1/3./sigma_t_./sigma_t_, ones(size(sigma_t_))) + alpha * h;

end