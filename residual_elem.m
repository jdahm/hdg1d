function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
	 residual_elem(Q, U, L, xg, lisb, risb, lbd, rbd, dx, td, fd, sd, qd)

  nnq = size(qd.qPhi, 2);
  nnv = size(qd.vPhi, 2);
  nnu = size(qd.uPhi, 2);
  nnw = size(qd.wPhi, 2);

  UH_Q = zeros([2,nnq]);
  UH_U = zeros([2,nnu]);
  UH_L = zeros([2,2]);

  if lisb
    q0 = qd.qPhi0*Q;
    u0 = qd.uPhi0*U;
    [ub, ub_q, ub_u] = boundary_state(q0, u0, fd, lbd);

    UH(1) = ub;
    UH_Q(1,:) = ub_q*qd.qPhi0;
    UH_U(1,:) = ub_u*qd.uPhi0;
  else
    UH(1,1) = L(1);
    UH_L(1,1) = 1.;
  end

  if risb
    q1 = qd.qPhi1*Q;
    u1 = qd.uPhi1*U;
    [ub, ub_q, ub_u] = boundary_state(q1, u1, fd, rbd);

    UH(2) = ub;
    UH_Q(2,:) = ub_q*qd.qPhi1;
    UH_U(2,:) = ub_u*qd.uPhi1;
  else
    UH(2,1) = L(2);
    UH_L(2,2) = 1.;
  end

  [Rq, Rq_Q, Rq_U, Rq_UH] = Rq_elem(Q, U, UH, dx, qd);
  Rq_Q = Rq_Q + Rq_UH*UH_Q;
  Rq_U = Rq_U + Rq_UH*UH_U;
  Rq_L = Rq_UH*UH_L;

  [Ru, Ru_Q, Ru_U, Ru_UH] = Ru_elem(Q, U, UH, xg, td, fd, sd, dx, lisb, risb, qd);
  Ru_Q = Ru_Q + Ru_UH*UH_Q;
  Ru_U = Ru_U + Ru_UH*UH_U;
  Ru_L = Ru_UH*UH_L;

  [Rl, Rl_Q, Rl_U, Rl_UH] = Rl_elem(Q, U, UH, fd, lisb, risb, qd);
  Rl_Q = Rl_Q + Rl_UH*UH_Q;
  Rl_U = Rl_U + Rl_UH*UH_U;
  Rl_L = Rl_UH*UH_L;
