function [dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  %
  % PURPOSE: Builds the HDG system K*L=R, solves it, then computes the updates to Q, U, L
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nelem*nn{q,u}]
  %   L : basis coefficients for trace [nelem-1]
  %   {l,r}bd : {left,right} boundary data [struct]
  %   md : mesh data [struct]
  %   td : time data [struct]
  %   fd : flux data [struct]
  %   sd : source data [struct]
  %   qd : quadrature data [struct]
  %
  % OUTPUTS:
  %   d{Q,U} : update for {Q,U} [nelem*nn{q,u}]
  %   dL : update for L [nelem-1]

  % 1. Build system
  [K, R] = hdg_system(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

  % 2. Solve reduced HDG system (linear solve)
  % dL = K\R;

  % 3. Post-linear solve
  % [dQ, dU] = hdg_post(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

function [K, R] = hdg_system(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [K, R] = hdg_system(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  %
  % PURPOSE: Loop over elements and add contribution to Jacobian
  %
  % NOTE: first and last elements need to be handled separately to
  % enforce the boundary condition
  %
  % INPUTS: (same as main function)
  %
  % OUTPUTS:
  %   K : HDG matrix [nelem-1,nelem-1]
  %   R : RHS vector [nelem-1]
  %

  % initialize R here, K will be built sparse at once at the end
  R = zeros([md.ne-1,1]);
  B = zeros([md.ne-1,3]);

  nnq = qd.pq+1;
  nnu = qd.pu+1;

  % first element
  QE = Q(1:nnq);
  UE = U(1:nnu);

  [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, lbd);
  UH(1) = ub;
  UH_QE(1) = ub_q;
  UH_UE(1) = ub_u;
  UH_LE(1) = 0.;

  if md.ne == 1
    [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, rbd);
    UH(2) = ub;
    UH_QE(2) = ub_q;
    UH_UE(2) = ub_u;
    UH_LE(2) = 0.;
    [~, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
    % no contribution to K
    % in this case, no global linear system is solved, only need post-linear solve routine
  else
    UH(2) = L(1);
    UH_QE(2) = 0.;
    UH_QE(2) = 0.;
    UH_LE(2) = 1.;
    [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
    B(1,2) = BK(2,2);
    R(1) = BR(2);
  end

  % interior elements
  for elem=2:md.ne-1
      QE = Q(1:nnq);
      UE = U(1:nnu);
      UH(1) = L(elem-1);
      UH_QE(1) = 0.;
      UH_UE(1) = 0.;
      UH_LE(1) = 1.;
      UH(2) = L(elem);
      UH_QE(2) = 0.;
      UH_UE(2) = 0.;
      UH_LE(2) = 1.;
      [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
      B(elem-1,2) = B(elem-1,2) + BK(1,1);
      B(elem-1,1) = B(elem-1,1) + BK(2,1);
      B(elem-1,3) = B(elem-1,3) + BK(1,2);
      B(elem,2)   = B(elem,2)   + BK(2,2);    
      R(elem-1:elem) = R(elem-1:elem) + BR;
  end

  % last element
  if md.ne > 1
    QE = Q(1:nnq);
    UE = U(1:nnu);

    UH(1) = L(nelem-1);
    UH_QE(1) = ub_q;
    UH_UE(1) = ub_u;
    UH_LE(1) = 0.;

    [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, rbd);
    UH(2) = ub;
    UH_QE(2) = ub_q;
    UH_UE(2) = ub_u;
    UH_LE(2) = 0.;

    [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
    B(md.ne-1,2) = B(md.ne-1,2) + BK(1,1);
    R(md.ne-1) = R(md.ne-1) + BR(1);
  end

  K = spdiags(B, [-1:1], md.ne-1, md.ne-1);

function [dQ, dU] = hdg_post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU] = hdg_post(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  %
  % PURPOSE: After linear solve for dL, loop over elements and obtain dQ and dU.
  %
  % NOTE: first and last elements need to be handled separately to
  % enforce the boundary condition
  %
  % INPUTS: (same as main function)
  %
  % OUTPUTS:
  %   d{Q,U} : update for {Q,U} [nelem*nn{q,u}]
  %

  nnq = qd.pq+1;
  nnu = qd.pu+1;

  % first element
  QE = Q(1:nnq);
  UE = U(1:nnu);

  [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, lbd);
  UH(1) = ub;
  UH_QE(1) = ub_q;
  UH_UE(1) = ub_u;
  UH_LE(1) = 0.;

  if nelem == 1
    [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, rbd);
    UH(2) = ub;
    UH_QE(2) = ub_q;
    UH_UE(2) = ub_u;
    UH_LE(2) = 0.;
    [~, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
    % no contribution to K
    % in this case, no global linear system is solved, only need post-linear solve routine
  else
    UH(2) = L(1);
    UH_QE(2) = 0.;
    UH_QE(2) = 0.;
    UH_LE(2) = 1.;
    [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
  end

  % interior elements
  for elem=2:md.ne-1
      QE = Q(1:nnq);
      UE = U(1:nnu);
      UH(1) = L(elem-1);
      UH_QE(1) = 0.;
      UH_UE(1) = 0.;
      UH_LE(1) = 1.;
      UH(2) = L(elem);
      UH_QE(2) = 0.;
      UH_UE(2) = 0.;
      UH_LE(2) = 1.;
      [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
  end

  % last element
  if md.ne > 1
    QE = Q(1:nnq);
    UE = U(1:nnu);

    UH(1) = L(nelem-1);
    UH_QE(1) = ub_q;
    UH_UE(1) = ub_u;
    UH_LE(1) = 0.;

    [ub, ub_q, ub_u] = boundary_state(QE, UE, td, fd, rbd);
    UH(2) = ub;
    UH_QE(2) = ub_q;
    UH_UE(2) = ub_u;
    UH_LE(2) = 0.;

    [BK, BR] = elem_contribution(QE, UE, UH, UH_QE, UH_UE, UH_LE, md.dx, td, fd, sd, qd);
  end

function [BK, BR] = elem_contribution(Q, U, UH, UH_q, UH_u, UH_l, dx, td, fd, sd, qd)

  % calculate residual and linearization blocks
  if fd.q_present
    [Rq, Rq_Q, Rq_U, Rq_UH] = Rq_elem(Q, U, UH, dx, qd);
  end
  size(Rq_UH)
  size(UH_Q)
  size(Rq_Q)
  Rq_Q = Rq_Q + bsxfun(@times, Rq_UH, UH_q);
  Rq_U = Rq_U + bsxfun(@times, Rq_UH, UH_q);
  Rq_L = bsxfun(@times, Rq_UH, UH_L);
  Rq

  [Ru, Ru_Q, Ru_U, Ru_UH] = Ru_elem(Q, U, UH, td, fd, sd, dx, qd);
  Ru_Q = Ru_Q + bsxfun(@times, Ru_UH, UH_Q);
  Ru_U = Ru_U + bsxfun(@times, Ru_UH, UH_U);
  Ru_L = bsxfun(@times, Ru_UH, UH_L);
  Ru

  [Rl, Rl_Q, Rl_U, Rl_UH] = Rl_elem(Q, U, UH, fd, qd);
  Rl_Q = Rl_Q + bsxfun(@times, Rl_UH, UH_Q);
  Rl_U = Rl_U + bsxfun(@times, Rl_UH, UH_U);
  Rl_L = bsxfun(@times, Rl_UH, UH_L);
  Rl
  BK = zeros([2,2]);
  BR = zeros([2,1]);
  return

  % initialize matrix and vector
  % BK = D
  % BR = H
  BK = Rl_L;
  BR = -Rl; % want -R on RHS

  % pre-process A matrix
  [A11, A12, A21, A22] = preprocess_A(Rq_Q, Rq_U, Ru_Q, Ru_U, fd.q_present);

  % subtract CA^{-1}B from BK
  % loop over faces (in 1D, left then right)
  [iABQ, iABU] = apply_Ainv(A11, A12, A21, A22, Rq_L, Ru_L);
  BK = BK - (Rl_Q*iABQ + Rl_U*iABU);

  % subtract C*A^{-1}*F from BR
  % this time A^{-1} multiplies a vector (F) from the left, so don't
  % need to loop over faces (columns) of F
  [iAFQ, iAFU] = apply_Ainv(A11, A12, A21, A22, Rq, Ru);
  BR = BR - (Rl_Q*iAFQ + Rl_U*iAFU);

function [A11, A12, A21, A22] = preprocess_A(RQQ, RQU, RUQ, RUU, q_present)

  A11 = RQQ;
  A12 = RQU;
  A21 = RUQ;
  A22 = RUU;

  if q_present % needed because A11 should not be invertible unless q is present
    A12 = A11\A12;
    A22 = A11 - A21*A12;
    A21 = A22\A21;
  end

function [iAQ, iAU] = apply_Ainv(A11, A12, A21, A22, Q, U)

  if q_present
    % needed again, A11 should not be invertible unless q is present
    BQ = A11\Q;
  end
  WA = -A22\U;
  if q_present
    % TODO could set Q and BQ to zero to have the same effect as this
    % if statement
    WA = WA + A21*Q;
    iAQ = BQ + A12*WA;
  end
  iAU = -WA;
