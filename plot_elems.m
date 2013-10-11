function h = plot_elems(xs, xe, xn, U, np)
  % function h = plot_elems(xs, xe, xn, U, np)
  %
  % PURPOSE: Plots a scalar quantity over the elements in the mesh
  %
  % INPUTS:
  %   xs : global position of left side of mesh
  %   xe : global position of right side of mesh
  %   xn : points for the interpolating polynomial basis
  %   U : array of polyomial coefficients
  %   np : number of plotting points per element
  %
  % OUTPUTS:
  %   h : plot handle
  %

  xp = linspace(0, 1, np);
  Phi = basis(xn, xp');

  nn = size(Phi, 2);

  ne = length(U)/nn;

  dx = (xe-xs)/ne;

  hold on;
  for elem=1:ne
      y = Phi*U((elem-1)*nn+1:elem*nn);
      x = xs+(elem-1)*dx+dx*xp;
      h = plot(x,y,'-');
  end
  hold off;
  
