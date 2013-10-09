function [fs, fs_q, fs_u, fs_uh] = flux_stab(q, u, uh, n, fd)
  % function [fs, fs_q, fs_u, fs_uh] = flux_stab(q, u, uh, n, fd)
  %
  % PURPOSE: Calculates the flux stabilization term.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [1]
  %   u : scalar [1]
  %   uh : trace value at interface [1]
  %   n : normal vector [1]
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   fs : stabilization [1]
  %   fs_{q,u,uh} : derivative wrt {q,u,uh}
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

  fs   = tau*(u-uh);
  fs_q = 0.;
  fs_u = tau*1.0;
  fs_uh = -tau*1.0;
