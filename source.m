function [s, s_q, s_u] = source(q, u, sd)
  % function [s, s_q, s_u] = source(q, u, sd)
  %
  % PURPOSE: Computes the source term at an arbitrary number of points.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [np]
  %   u : scalar [np]
  %   sd : source data [struct]
  %
  % OUTPUTS:
  %   s : source [np]
  %   s_q : derivative of function wrt q [np]
  %   s_u : derivative of function wrt u [np]
  %

  s = sd.a*u + sd.b*q;
  s_q = sd.b*ones(size(q));
  s_u = sd.a*ones(size(u));
