function [s, s_q, s_u] = source(q, u, xg, sd)
  % function [s, s_q, s_u] = source(q, u, xg, sd)
  %
  % PURPOSE: Computes the source term at an arbitrary number of points.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [np]
  %   u : scalar [np]
  %   xglob : global positions of quadrature nodes [np]
  %   sd : source data [struct]
  %
  % OUTPUTS:
  %   s : source [np]
  %   s_q : derivative of function wrt q [np]
  %   s_u : derivative of function wrt u [np]
  %

  if sd.present
    switch (sd.type)
      case 'ms_quadratic'
	% u=2x-x^2 on [0,2]
	s = 2-2*xg+2;
	s_q = zeros(size(q));
	s_u = zeros(size(u));
      otherwise
	error('unknown source term');
    end
  else
      s = zeros(size(u));
      s_q = zeros(size(q));
      s_u = zeros(size(u));
  end
