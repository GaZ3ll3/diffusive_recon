global p

p = radiative();

p.sigma_a_ = p.sigma_a_0;
for i = 1 : 10
    p.LBFGS('obj', 'grad', 'callback_radiative', p.sigma_a_);
%     p.optimize(p.sigma_a_);
end
p.plot();
