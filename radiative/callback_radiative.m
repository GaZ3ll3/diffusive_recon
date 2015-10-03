function  callback_radiative(t, f, x)
global p
fprintf('%3d\t\t%0.8g\t\t%0.8g \n', t, f, norm(x - p.sigma_a)/norm(p.sigma_a));
% numofnodes = size(p.fem.Promoted.nodes, 2);
% trisurf(p.fem.TriMesh', p.fem.Promoted.nodes(1,1:numofnodes), ...
% p.fem.Promoted.nodes(2, 1:numofnodes), x(1:numofnodes),'EdgeColor', 'None');shading interp;
% colorbar;colormap jet;
% pause(0.2);
end

