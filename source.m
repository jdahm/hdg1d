function [s, s_q, s_u] = source(q, u, xg, fd, sd)
  % function [s, s_q, s_u] = source(q, u, xg, fd, sd)
  %
  % PURPOSE: Computes the source term at an arbitrary number of points.
  %
  % INPUTS:
  %   q : gradient of scalar (derivative in 1D: u_x) [np]
  %   u : scalar [np]
  %   xglob : global positions of quadrature nodes [np]
  %   fd : flux data [struct]
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
	% S = -a*u_x + b*u_xx = -a*(2-2*x)+b*(-2) = -2*(a+b)+2*x*a
	s = -2*(fd.a+fd.b)+2*fd.a*xg;
	s_q = zeros(size(q));
	s_u = zeros(size(u));
      case 'linear'
	% (b=0) s=sd.c*u -> u=exp(x/a)
	s = sd.a*u + sd.b*q;
	s_q = sd.b*ones(size(q));
	s_u = sd.a*ones(size(u));
      otherwise
	error('unknown source term');
    end
  else
      s = zeros(size(u));
      s_q = zeros(size(q));
      s_u = zeros(size(u));
  end
