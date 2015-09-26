function [r, g] = reg2(sigma_a_)

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

