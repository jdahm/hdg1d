function [fs, fs_u, fs_uh] = flux_stab(u, uh, n, fd)
  % function [fs, fs_u, fs_uh] = flux_stab(u, uh, n, fd)
  %
  % PURPOSE: Calculates the flux stabilization
  %
  % NOTE: term is NOT dotted with normal here
  %
  % INPUTS:
  %   u : scalar [1]
  %   uh : trace value at interface [1]
  %   n : normal vector [1]
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   fs : stabilization [1]
  %   fs_{u,uh} : derivative wrt {u,uh}
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
  tau = tau*fd.c;

  fs   = tau*(u-uh);
  fs_u = tau*1.0;
  fs_uh = -tau*1.0;
