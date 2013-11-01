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

  [tau, tau_u, tau_uh] = stab(u, uh, n, fd);

  fs   = tau*(u-uh);
  fs_u = tau + tau_u*u;
  fs_uh = -(tau + tau_uh*uh);
