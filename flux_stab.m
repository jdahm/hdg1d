function [fs, fs_q, fs_u, fs_uh] = flux_stab(q, u, uh, fd)
  % function [fs, fs_q, fs_u, fs_uh] = flux_stab(q, u, uh, fd)
  %
  % PURPOSE: Calculates the flux stabilization term.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [1]
  %   u : scalar [1]
  %   uh : trace value at interface [1]
  %   fd : flux data [struct]
  %
  % OUTPUTS:
  %   fs : stabilization [1]
  %   fs_{q,u,uh} : derivative wrt {q,u,uh}
  %

  tau = abs(fd.a);

  if fd.q_present
     tau = tau + abs(fd.b)/fd.vl;
  end

  fs   = tau*(u-uh);
  fs_q = zeros(size(q));
  fs_u = tau*ones(size(u));
  fs_uh = -tau*ones(size(uh));
