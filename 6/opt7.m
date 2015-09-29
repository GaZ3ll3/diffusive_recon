global p

p = opt6();

p.LBFGS('obj', 'grad', 'callback6');

p.plot();

