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
  dL = K\R;

  % 3. Post-linear solve
  %[dQ, dU] = hdg_post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd);
  dQ = zeros(size(Q));
  dU = zeros(size(U));

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
  LE(1,1) = 0.;

  lisb = true;
  if md.ne ~= 1
    risb = false;
    LE(2,1) = L(1);
  else
    risb = true;
    LE(2,1) = 0.;
  end

  [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_belem(QE, UE, LE, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  if md.ne ~= 1
    B(1,2) = BK(2,2);
    R(1) = BR(2);
  end

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_ielem(QE, UE, LE, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    B(elem-1,2) = B(elem-1,2) + BK(1,1);
    B(elem-1,1) = B(elem-1,1) + BK(2,1);
    B(elem-1,3) = B(elem-1,3) + BK(1,2);
    B(elem,2)   = B(elem,2)   + BK(2,2);    
    R(elem-1:elem) = R(elem-1:elem) + BR;
  end

  % last element
  if md.ne > 1
    QE = Q((md.ne-1)*nnq+1:md.ne*nnq);
    UE = U((md.ne-1)*nnu+1:md.ne*nnu);
    LE = [L(md.ne-1); 0.];
    LE
    lisb = false;
    risb = true;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_belem(QE, UE, LE, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    B(md.ne-1,2) = B(md.ne-1,2) + BK(1,1);
    R(md.ne-1) = R(md.ne-1) + BR(1);
  end

  K = spdiags(B, -1:1, md.ne-1, md.ne-1);


function [dQ, dU] = hdg_post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU] = hdg_post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  %
  % PURPOSE: After linear solve for dL, loop over elements and obtain dQ and dU.
  %
  % NOTE: first and last elements need to be handled separately to
  % enforce the boundary condition
  %
  % INPUTS: (same as main function)
  %   dL : Update to L [nelem-1]
  %
  % OUTPUTS:
  %   d{Q,U} : update for {Q,U} [nelem*nn{q,u}]
  %

  dQ = zeros(size(Q));
  dU = zeros(size(U));

  nnq = qd.pq+1;
  nnu = qd.pu+1;

  % first element
  QE = Q(1:nnq);
  UE = U(1:nnu);
  LE(1,1) = 0.;
  dLE(1,1) = 0.;

  lisb = true;
  if md.ne ~= 1
    risb = false;
    LE(2,1) = L(1);
    dLE(2,1) = dL(1);
  else
    risb = true;
    LE(2,1) = 0.;
    dLE(2,1) = 0.;
  end

  [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_belem(QE, UE, LE, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  dQ(1:nnq) = dQE;
  dU(1:nnu) = dUE;

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);
    dLE = dL(elem-1:elem);

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_ielem(QE, UE, LE, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((elem-1)*nnq+1:elem*nnq) = dQE;
    dU((elem-1)*nnu+1:elem*nnu) = dUE;
  end

  % last element
  if md.ne > 1
    QE = Q((md.ne-1)*nnq+1:md.ne*nnq);
    UE = U((md.ne-1)*nnu+1:md.ne*nnu);
    LE  = [L(md.ne-1); 0.];
    dLE = [L(md.ne-1); 0.];

    lisb = false;
    risb = true;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_belem(QE, UE, LE, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((md.ne-1)*nnq+1:md.ne*nnq) = dQE;
    dU((md.ne-1)*nnu+1:md.ne*nnu) = dUE;
  end


function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_belem(Q, U, L, lisb, risb, lbd, rbd, dx, td, fd, sd, qd)

  nnu = size(qd.qPhi, 2);
  nnq = size(qd.uPhi, 2);

  UH_Q = zeros([2,nnq]);
  UH_U = zeros([2,nnu]);
  UH_L = zeros([2,2]);

  if lisb
    q0 = qd.qPhi0*Q;
    u0 = qd.uPhi0*U;
    [ub, ub_q, ub_u] = boundary_state(q0, u0, td, fd, lbd);

    UH(1) = ub;
    UH_Q(1,:) = ub_q*qd.qPhi0;
    UH_U(1,:) = ub_u*qd.uPhi0;
    UH_L(1,1) = 0.;
  else
    UH(1) = L(1);
    UH_Q(1,:) = zeros(size(qd.qPhi1));
    UH_U(1,:) = zeros(size(qd.uPhi1));
    UH_L(1,1) = 1.;
  end

  if risb
    q1 = qd.qPhi1*Q;
    u1 = qd.uPhi1*U;
    [ub, ub_q, ub_u] = boundary_state(q1, u1, td, fd, rbd);

    UH(2) = ub;
    UH_Q(2,:) = ub_q*qd.qPhi1;
    UH_U(2,:) = ub_u*qd.uPhi1;
    UH_L(2,2) = 0.;
  else
    UH(2) = L(2);
    UH_Q(2,:) = zeros(size(qd.qPhi1));
    UH_U(2,:) = zeros(size(qd.uPhi1));
    UH_L(2,2) = 1.;
  end

  [Rq, Rq_Q, Rq_U, Rq_UH] = Rq_elem(Q, U, UH, dx, qd);
  Rq_Q = Rq_Q + Rq_UH*UH_Q;
  Rq_U = Rq_U + Rq_UH*UH_U;
  Rq_L = Rq_UH*UH_L;

  [Ru, Ru_Q, Ru_U, Ru_UH] = Ru_elem(Q, U, UH, td, fd, sd, dx, qd);
  Ru_Q = Ru_Q + Ru_UH*UH_Q;
  Ru_U = Ru_U + Ru_UH*UH_U;
  Ru_L = Ru_UH*UH_L;

  [Rl, Rl_Q, Rl_U, Rl_UH] = Rl_elem(Q, U, UH, fd, lisb, risb, qd);
  Rl_Q = Rl_Q + Rl_UH*UH_Q;
  Rl_U = Rl_U + Rl_UH*UH_U;
  Rl_L = Rl_UH*UH_L;

function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = residual_ielem(Q, U, L, dx, td, fd, sd, qd)

  [Rq, Rq_Q, Rq_U, Rq_L] = Rq_elem(Q, U, L, dx, qd);

  [Ru, Ru_Q, Ru_U, Ru_L] = Ru_elem(Q, U, L, td, fd, sd, dx, qd);

  [Rl, Rl_Q, Rl_U, Rl_L] = Rl_elem(Q, U, L, fd, false, false, qd);


function [dQ, dU] = qu_backsolve(dL, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % pre-process A matrix
  [AQQ, AQU, AUQ, AUU] = preprocess_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  Rq = -Rq - Rq_L*dL;
  Ru = -Ru - Ru_L*dL;
  [dQ, dU] = apply_Ainv(AQQ, AQU, AUQ, AUU, Rq, Ru, q_present);


function [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % initialize matrix and vector
  % BK = D
  % BR = H
  BK = Rl_L;
  BR = -Rl; % want -R on RHS

  % pre-process A matrix
  [AQQ, AQU, AUQ, AUU] = preprocess_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  % subtract CA^{-1}B from BK
  % loop over faces (in 1D, left then right)
  [iABQ, iABU] = apply_Ainv(AQQ, AQU, AUQ, AUU, Rq_L, Ru_L, q_present);
  BK = BK - (Rl_Q*iABQ + Rl_U*iABU);

  % subtract C*A^{-1}*F from BR, where F=-[Rq;Ru]
  % this time A^{-1} multiplies a vector (F) from the left, so don't
  % need to loop over faces (columns) of F
  [iAFQ, iAFU] = apply_Ainv(AQQ, AQU, AUQ, AUU, -Rq, -Ru, q_present);
  BR = BR - (Rl_Q*iAFQ + Rl_U*iAFU);

function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)

  if q_present % needed because AQQ should not be invertible unless q is present
    % AQU = AQQ^{-1}*AQU
    AQU = AQQ\AQU;
    % AUU -= AUQ*AQU
    AUU = AUU-AUQ*AQU;
    % AUQ = AUU^{-1}*AUQ
    AUQ = AUU\AUQ;
  end


function [Q, U] = apply_Ainv(AQQ, AQU, AUQ, AUU, Q, U, q_present)

  if q_present
    % needed again, AQQ should not be invertible unless q is present
     % Q = AQQ^{-1}*Q
    Q = AQQ\Q;
  end
  % W = -AUU^{-1}*U
  W = -AUU\U;
  if q_present
    % TODO could set Q and BQ to zero to have the same effect as this
    % if statement
    % W += AUQ*Q
    W = W+AUQ*Q;
    % Q += AQU*W
    Q = Q+AQU*W;
  end
  % U = -W
  U = -W;
