function [tau, tau_u, tau_uh] = stab(u, uh, n, fd)
  % function [tau, tau_u, tau_uh] = stab(u, uh, n, fd)
  %
  % PURPOSE: Calculates the stabilization $tau$
  %
  % INPUTS:
  %   u : element solution
  %   uh : trace of u
  %   n : normal vector (1 for right face, -1 for left face)
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   tau : stabilization term
  %   tau_{u,uh} : linearization of tau wrt. {u,trace}
  %

  switch fd.stab_type
    case 'centered'
      tau = abs(fd.a);
      if fd.q_present
	tau = tau + fd.b/fd.vl;
      end
    case 'upwind'
      tau = 0.5*(fd.a*n+abs(fd.a*n));
      if fd.q_present
	tau = tau + fd.b/fd.vl;
      end
  end
  tau = fd.c*tau;
  tau_u = zeros(size(u));
  tau_uh = zeros(size(uh));
