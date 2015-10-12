global p

p = radiative3();

p.LBFGS('obj3', 'grad3', 'callback_radiative3');

p.plot();
