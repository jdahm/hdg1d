function [PQ, PU] = postprocess_dirichlet(Q, U, L, xnq, pxnq, xnu, pxnu, lbd, rbd, md, td, fd, sd)

  nnq = length(xnq);
  pnnq = length(pxnq);
  nnu = length(xnu);
  pnnu = length(pxnu);

  pq = nnq-1;
  pu = nnu-1;

  qd = quad_data(pxnq, pxnq, pxnu, pxnu);

  PQ = zeros([md.ne*pnnq,1]);
  PU = zeros([md.ne*pnnu,1]);

  TQ = basis(xnq, pxnq);
  TU = basis(xnu, pxnu);

  for elem=1:md.ne
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);

    % inject q and u to the post-processed space
    PQE = TQ*QE;
    PUE = TU*UE;

    if elem ~= 1
      uh(1) = L(elem-1);
    else
      q0 = qd.qPhi0*PQE;
      u0 = qd.uPhi0*PUE;
      [uh(1), ~, ~] = boundary_state(q0, u0, fd, lbd);
    end
    if elem ~= md.ne
      uh(2) = L(elem);
    else
      q1 = qd.qPhi1*PQE;
      u1 = qd.uPhi1*PUE;
      [uh(2), ~, ~] = boundary_state(q1, u1, fd, rbd);
    end

    xg = (elem-1)*md.dx+qd.x*md.dx;

    [PQE, PUE] = postprocess_elem(PQE, PUE, uh, xg, md.dx, td, fd, sd, qd);

    PQ((elem-1)*pnnq+1:elem*pnnq) = PQE;
    PU((elem-1)*pnnu+1:elem*pnnu) = PUE;

  end


function [Q, U] = postprocess_elem(Q, U, uh, xg, dx, td, fd, sd, qd)

  nnq = size(qd.qPhi, 2);
  nnu = size(qd.uPhi, 2);

  [Rq, Ru, ~, Rq_Q, Rq_U, ~, Ru_Q, Ru_U, ~, ~, ~, ~] = ...
  residual_elem(Q, U, uh, xg, false, false, 0, 0, dx, td, fd, sd, qd);

  dQU = [Rq_Q, Rq_U; Ru_Q, Ru_U] \ [-Rq; -Ru];
  Q = Q + dQU(1:nnq);
  U = U + dQU(nnq+1:nnq+nnu);
