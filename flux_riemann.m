function [f, f_q, f_u, f_uh] = flux_riemann(q, u, uh, n, fd)
  % function [f, f_q, f_u, f_uh] = flux_bc(q, u, uh, n, fd)
  %
  % PURPOSE: Computes the flux dotted with normal.
  % This version uses a Riemann solver to properly 
  % upwind the advective part.
  %
  % INPUTS:
  %   {q,u,uh} : {grad(u),u,boundary value} at interface
  %   n : normal vector
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   f : F.n
  %   f_{q,u,uh} : linearization of flux wrt {q,u,uh}
  %

  % simple in this case:
  % if a*n > 0 then the interface state ("*-state") is given by the
  % left state, otherwise it's the right state
  if fd.a*n > 0
    [f, f_q, f_u] = flux(q, u, fd);
    f_uh = zeros(size(uh));
  else
    [f, f_q, f_uh] = flux(q, uh, fd);
    f_u = zeros(size(u));
  end

  f = n*f;
  f_q = n*f_q;
  f_u = n*f_u;
