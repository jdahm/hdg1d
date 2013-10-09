function [f, f_q, f_u] = flux(q, u, fd)
  % function [f, f_q, f_u] = flux(q, u, fd)
  %
  % PURPOSE: Computes the flux and derivatives at an arbitrary number of points.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [np]
  %   u : scalar [np]
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   f : flux [np]
  %   f_q : derivative of function wrt q
  %   f_u : derivative of function wrt u
  %

  f   = fd.a*u - fd.b*q;
  f_u = fd.a*ones(size(u));
  f_q = -fd.b*ones(size(q));
