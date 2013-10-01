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

  md.ne = 1;
  md.xs = 0.;
  md.xe = 1.;
  md.dx = (xe-xs)/ne;

