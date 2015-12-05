function [f] = data2(sigma_a_)

global fem u load sigma_s sigma_a Q counter sigma_rate M alpha

sigma_t_ = sigma_a_ + sigma_s;

T = fem.assems((1/3)./sigma_t_) + fem.assema(sigma_a_) + 0.5 * Q;

v = T\load;

[r, ~] = reg2(sigma_a_);

% f = 0.5 * (u - v)' * M *  (u - v) + alpha * r;
f = v' * (log(v) - log(u)) + sum(u - v) + alpha * r;

sigma_rate(counter) = norm(sigma_a_ - sigma_a)/norm(sigma_a);
counter = counter + 1;


end