function [f, f_q, f_u, f_uh] = flux_oneside(q, u, uh, fd)
  % function [f, f_q, f_u, f_uh] = flux_oneside(q, u, uh, fd)
  %
  % PURPOSE: Computes the one-sided flux at a single point.
  %
  % NOTE: This function assumes that fd.stab is independent of both u
  % and the trace uh.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [1]
  %   u : scalar [1]
  %   uh : trace value at interface [1]
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   f : flux [1]
  %   f_q : derivative of function wrt q
  %   f_u : derivative of function wrt u
  %   f_uh : derivative of function wrt uh
  %

  [h, h_q, h_u] = flux(q,u,fd);
  f = h + fd.stab*(u-uh);
  f_q = h_q;
  f_u = h_u + fd.stab;
  f_uh = -1.0*fd.stab;
