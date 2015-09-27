function callback (t, f, x)
global sigma_a
  fprintf('%5d\t\t%0.8g\t\t%6.8f\n', t, f, norm(x - sigma_a));
end

