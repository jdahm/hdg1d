function [Q, U, L] = initialize(pq, pu, nelem)
  % function [Q, U, L] = initialize(pq, pu)
  %
  % PURPOSE: A shortcut to set the solution to zeros.
  %
  % INPUTS:
  %   pu : order of the polynomial for the scalar on each element
  %   pq : order of the polynomial for the gradient of the scalar on each element
  %   nelem : number of elements
  %
  % OUTPUTS:
  %   {Q,U,L} : basis coefficients
  %   

  nnq = pq+1;
  nnu = pu+1;

  Q = zeros([nnq*nelem,1]);

  U = zeros([nnu*nelem,1]);

  % no trace on boundary faces
  L = zeros([nelem-1,1]);
