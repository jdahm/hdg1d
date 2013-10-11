function h = plot_traces(xs, xe, L)
  % function h = plot_traces(xs, xe, L)
  %
  % PURPOSE: Plots traces on interfaces of the mesh
  %
  % INPUTS:
  %   xs : global position of left side of mesh
  %   xe : global position of right side of mesh
  %   L : trace values
  %
  % OUTPUTS:
  %   h : plot handle
  %

  ne = length(L)+1;

  dx = (xe-xs)/ne;

  hold on;
  for elem=1:ne-1
      x = elem*dx;
      y = L(elem);
      h = plot(x, y, '.');
  end
  hold off;
