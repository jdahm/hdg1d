function md = mesh(xs, xe, ne)
  % function md = mesh(xs, xe, ne)
  %
  % INPUTS:
  %   xs : global position of left side of mesh
  %   xe : global position of right side of mesh
  %   ne : number of mesh elements
  %
  % OUTPUTS:
  %   md : mesh data
  %

  md.ne = ne;
  md.xs = xs;
  md.xe = xe;
  md.dx = (xe-xs)/ne;

