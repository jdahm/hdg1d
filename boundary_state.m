function [ub, ub_q, ub_u] = boundary_state(q, u, td, fd, bd)
  % function [ub, ub_q, ub_u] = boundary_state(q, u, td, fd, bd)
  %
  % PURPOSE: Computes the boundary state ub(q,u) and linearization
  % at an arbitrary number of points.
  %
  % INPUTS:
  %   q : grad(u) (= u_x for 1D) at points [np]
  %   u : scalar at points [np]
  %   td : time data [struct]
  %   fd : flux data [struct]
  %   bd : boundary data [struct]
  %
  % OUTPUTS:
  %   ub : boundary state [np]
  %   ub_q : linearization of boundary state wrt q [np]
  %   ub_u : linearization of boundary state wrt u [np]
  %

  switch bd.type
    case 'd'
      ub = bd.data(1);
      ub_q = 0.;
      ub_u = 0.;
    otherwise
      error('unknown boundary type');
  end
