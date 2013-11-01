function [QP, UP] = postprocess_neumann(Q, U, L, qbasis, pq, ubasis, pu, lbd, rbd, md, td, fd, sd)

  xnq = create_nodes(pq, qbasis);
  xnq1 = create_nodes(pq+1, qbasis);
  xnqm1 = create_nodes(pq-1, qbasis);

  xnu = create_nodes(pu, ubasis);
  xnu1 = create_nodes(pu+1, ubasis);

  qorder = max([pq+1,pu+1])+1;

  pqd = struct;

  [pqd.x, pqd.w] = lgwt(qorder, 0, 1);

  pqd.qPhiqm1 = basis(xnqm1, pqd.x);
  pqd.qPhiq = basis(xnq, pqd.x);
  pqd.qPhiq_0 = basis(xnq, 0.);
  pqd.qPhiq_1 = basis(xnq, 1.);

  pqd.qPhiu1 = basis(xnu1, pqd.x);
  pqd.qPhiq1 = basis(xnq1, pqd.x);
  pqd.qPhiq1_0 = basis(xnq1, 0.);
  pqd.qPhiq1_1 = basis(xnq1, 1.);

  pqd.qPhiu = basis(xnu, pqd.x);
  pqd.qgPhiu1 = gbasis(xnu1, pqd.x);

  pqd.uPhiu = basis(xnu, pqd.x);
  pqd.uPhiu1 = basis(xnu1, pqd.x);

  pqd.uPhiu_0 = basis(xnu, 0.);
  pqd.uPhiu_1 = basis(xnu, 1.);

  pqd.ugPhiu1 = gbasis(xnu1, pqd.x);

  nnq = length(xnq);
  nnu = length(xnu);
  nnq1 = length(xnq1);
  nnu1 = length(xnu1);

  QP = zeros([md.ne*nnq1,1]);
  UP = zeros([md.ne*nnu1,1]);

  xplot = linspace(0, 1, 100)';
  Phiqplot = basis(xnq, xplot);
  Phiuplot = basis(xnu, xplot);
  Phiq1plot = basis(xnq1, xplot);
  Phiu1plot = basis(xnu1, xplot);

  figure(99); clf; hold on;
  figure(100); clf; hold on;

  for elem=1:md.ne
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    
    if elem ~= 1
      uh(1) = L(elem-1);
    else
      q0 = pqd.qPhiq_0*QE;
      u0 = pqd.uPhiu_0*UE;
      [uh(1), ~, ~] = boundary_state(q0, u0, fd, lbd);
    end
    if elem ~= md.ne
      uh(2) = L(elem);
    else
      q1 = pqd.qPhiq_1*QE;
      u1 = pqd.uPhiu_1*UE;
      [uh(2), ~, ~] = boundary_state(q1, u1, fd, rbd);
    end

    [QPE, UPE] = postprocess_elem(QE, UE, uh, fd, pqd, md.dx);
    xg = (elem-1)*md.dx+xplot*md.dx;

    figure(99);
    plot(xg, Phiqplot*QE, '-');
    plot(xg, Phiq1plot*QPE, '--');
    figure(100);
    plot(xg, Phiuplot*UE, '-');
    plot(xg, Phiu1plot*UPE, '--');

    QP((elem-1)*nnq1+1:elem*nnq1) = QPE;
    UP((elem-1)*nnu1+1:elem*nnu1) = UPE;

  end

  figure(99); hold off;
  figure(100); hold off;


function [QP, UP] = postprocess_elem(Q, U, uh, fd, pqd, dx)

  n0 = -1.;
  n1 = 1.;

  dw = diag(pqd.w);

  A = dx*pqd.qPhiqm1'*dw*pqd.qPhiq1;
  A0 = pqd.qPhiq1_0*n0;
  A1 = pqd.qPhiq1_1*n1;

  M = [A; A0; A1];

  b = dx*pqd.qPhiqm1'*dw*(pqd.qPhiq*Q);

  q0 = pqd.qPhiq_0*Q;
  u0 = pqd.uPhiu_0*U;

  q1 = pqd.qPhiq_1*Q;
  u1 = pqd.uPhiu_1*U;

  [tau0, ~, ~] = stab(u0, uh(1), n0, fd);
  [tau1, ~, ~] = stab(u1, uh(2), n1, fd);

  b0 = q0*n0 - tau0*(u0-uh(1))/fd.b;
  b1 = q1*n1 - tau1*(u1-uh(2))/fd.b;

  r = [b; b0; b1];

  QP = M\r;

  AP = dx*(pqd.ugPhiu1/dx)'*dw*pqd.ugPhiu1/dx;
  bP = dx*(pqd.ugPhiu1/dx)'*dw*pqd.qPhiq1*QP;

  cAP = dx*pqd.w'*pqd.uPhiu1;
  cbP = dx*pqd.w'*(pqd.uPhiu*U);

  UP = [AP; cAP]\[bP; cbP];
