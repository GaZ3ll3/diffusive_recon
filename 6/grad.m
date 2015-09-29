function g = grad(x)

global p

[~, g] = p.objective_gradient(x);

end

