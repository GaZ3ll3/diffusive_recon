function [r, g] = reg3(sigma_a_)

global fem

m =  fem.assems(1);
r = 0.5 * sigma_a_' * m * sigma_a_;
g = m * sigma_a_;


end