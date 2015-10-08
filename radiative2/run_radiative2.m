global p

p = radiative2();

p.LBFGS('obj2', 'grad2', 'callback_radiative2');

p.plot();
