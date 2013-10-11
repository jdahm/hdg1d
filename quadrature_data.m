function qd = quadrature_data(xnq, xnu, xnv, xnw)
  % function qd = quadrature_data(xnq, xnu, xnv, xnw)
  %
  % PURPOSE: Creates quadrature data struct used for integration of the bilinear forms
  %
  % INPUTS:
  %   xn{q,u,v,w} : nodal positions of roots for the interpolating polynomial basis defining {q,u,v} [nn{q,u,v,w}]
  %
  % OUTPUTS:
  %   qd : quadrature data struct  (see source code for members)
  %

  % create struct
  qd = struct;
  
  % store nodes
  qd.xnq = xnq;
  qd.xnv = xnv;
  qd.xnu = xnu;
  qd.xnw = xnw;

  % polynomial orders
  qd.pq = length(xnq)-1;
  qd.pv = length(xnv)-1;

  qd.pu = length(xnu)-1;
  qd.pw = length(xnw)-1;

  % quadrature order sufficient to integrate 2*p
  qd.order = max([qd.pq, qd.pv, qd.pu, qd.pw])+1;

  % quadrature points and weights
  [qd.x, qd.w] = lgwt(qd.order, 0, 1);
  qd.dw = diag(qd.w);

  % bases
  qd.vPhi  = basis(xnv, qd.x);
  qd.vPhi0 = basis(xnv, 0.0);
  qd.vPhi1 = basis(xnv, 1.0);

  qd.wPhi  = basis(xnw, qd.x);
  qd.wPhi0 = basis(xnw, 0.0);
  qd.wPhi1 = basis(xnw, 1.0);

  qd.qPhi  = basis(xnq, qd.x);
  qd.qPhi0 = basis(xnq, 0.0);
  qd.qPhi1 = basis(xnq, 1.0);

  qd.uPhi  = basis(xnu, qd.x);
  qd.uPhi0 = basis(xnu, 0.0);
  qd.uPhi1 = basis(xnu, 1.0);

  % gradients wrt reference space
  qd.vGPhi = gbasis(xnv, qd.x);
  qd.wGPhi = gbasis(xnw, qd.x);
