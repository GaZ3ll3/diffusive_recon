function f = obj(x)

global p

[f, ~] = p.objective_gradient(x);

end


