function [dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  %
  % PURPOSE: Builds the HDG system K*L=F, solves it, then computes the updates to Q, U, L
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

  % TEMP ping test
  % pass = hdg_ping(Q, U, L, lbd, rbd, md, td, fd, sd, qd)

  % 1. Build system (K*L=F) (F=-Residual)
  [K, F] = hdg_system(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

  % 2. Solve reduced HDG system (linear solve)
  dL = K\F;

  % 3. Post-linear solve
  [dQ, dU] = hdg_post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd);


% -----------------------------------------------------------------------
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
  
  xglob = qd.x*md.dx;

  [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
  hdg_residual_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				 Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  if md.ne ~= 1
    B(1,2) = BK(2,2);
    R(1) = BR(2);
  end

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);

    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_ielem(QE, UE, LE, xg, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				   Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    % self-block
    B(elem-1,2) = B(elem-1,2) + BK(1,1);
    B(elem,2)   = B(elem,2)   + BK(2,2);
    % off-diagonal entries
    B(elem-1,1) = B(elem-1,1) + BK(2,1);
    B(elem,3) = B(elem,3) + BK(1,2);
    R(elem-1:elem) = R(elem-1:elem) + BR;
  end

  % last element
  if md.ne > 1
    QE = Q((md.ne-1)*nnq+1:md.ne*nnq);
    UE = U((md.ne-1)*nnu+1:md.ne*nnu);
    LE = [L(md.ne-1); 0.];
    lisb = false;
    risb = true;

    xg = (md.ne-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				   Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    B(md.ne-1,2) = B(md.ne-1,2) + BK(1,1);
    R(md.ne-1) = R(md.ne-1) + BR(1);
  end

  K = spdiags(B, -1:1, md.ne-1, md.ne-1);


% -----------------------------------------------------------------------
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

  xg = qd.x*md.dx;

  [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
  hdg_residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			    Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  dQ(1:nnq) = dQE;
  dU(1:nnu) = dUE;

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);
    dLE = dL(elem-1:elem);

    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_ielem(QE, UE, LE, xg, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			      Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((elem-1)*nnq+1:elem*nnq) = dQE;
    dU((elem-1)*nnu+1:elem*nnu) = dUE;
  end

  % last element
  if md.ne > 1
    QE = Q((md.ne-1)*nnq+1:md.ne*nnq);
    UE = U((md.ne-1)*nnu+1:md.ne*nnu);
    LE  = [L(md.ne-1); 0.];
    dLE = [dL(md.ne-1); 0.];

    lisb = false;
    risb = true;

    xg = (md.ne-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			      Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((md.ne-1)*nnq+1:md.ne*nnq) = dQE;
    dU((md.ne-1)*nnu+1:md.ne*nnu) = dUE;
  end


% -----------------------------------------------------------------------
function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
	 hdg_residual_elem(Q, U, L, xg, lisb, risb, lbd, rbd, dx, td, fd, sd, qd)

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


% -----------------------------------------------------------------------
function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
	 hdg_residual_ielem(Q, U, L, xg, dx, td, fd, sd, qd)
  % function [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
  %          hdg_residual_ielem(Q, U, L, xg, dx, td, fd, sd, qd)
  %
  % Companion function to hdg_residual_elem (see function for documentation)
  %

  [Rq, Rq_Q, Rq_U, Rq_L] = Rq_elem(Q, U, L, dx, qd);
  [Ru, Ru_Q, Ru_U, Ru_L] = Ru_elem(Q, U, L, xg, td, fd, sd, dx, false, false, qd);
  [Rl, Rl_Q, Rl_U, Rl_L] = Rl_elem(Q, U, L, fd, false, false, qd);


% -----------------------------------------------------------------------
function [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
					Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % initialize matrix and vector
  % BK = D
  % BR = H
  BK = Rl_L;
  BR = -Rl;

  % pre-process A matrix
  [AQQ, AQU, AUQ, AUU] = preprocess_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  % subtract C(A^{-1}B) from BK
  [iABQ, iABU] = apply_Ainv(AQQ, AQU, AUQ, AUU, Rq_L, Ru_L, q_present);
  BK = BK - Rl_U*iABU;
  if q_present
    BK = BK - Rl_Q*iABQ;
  end

  % subtract C(A^{-1}F) from BR, where F=-[Rq;Ru]
  [iAFQ, iAFU] = apply_Ainv(AQQ, AQU, AUQ, AUU, -Rq, -Ru, q_present);
  BR = BR - Rl_U*iAFU;
  if q_present
    BR = BR - Rl_Q*iAFQ;
  end

% -----------------------------------------------------------------------
function [dQ, dU] = qu_backsolve(dL, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				 Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % pre-process A matrix
  [AQQ, AQU, AUQ, AUU] = preprocess_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  % F - B*dL
  Rq = -Rq - Rq_L*dL;
  Ru = -Ru - Ru_L*dL;

  [dQ, dU] = apply_Ainv(AQQ, AQU, AUQ, AUU, Rq, Ru, q_present);


% -----------------------------------------------------------------------
function pass = hdg_ping(Q, U, L, lbd, rbd, md, td, fd, sd, qd)

  % assume passed
  pass = true;

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

  xglob = qd.x*md.dx;

  pass_elem = ping_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
  if ~pass_elem
     pass = false;
     return;
  end

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);

    xg = (md.ne-1)*md.dx+qd.x*md.dx;
    lisb = false;
    risb = false;

    pass_elem = ping_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
    if ~pass_elem
      pass = false;
      return;
    end

  end

  % last element
  if md.ne > 1
    QE = Q((md.ne-1)*nnq+1:md.ne*nnq);
    UE = U((md.ne-1)*nnu+1:md.ne*nnu);
    LE = [L(md.ne-1); 0.];
    lisb = false;
    risb = true;

    xg = (md.ne-1)*md.dx+qd.x*md.dx;

    pass_elem = ping_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
    if ~pass_elem
      pass = false;
      return;
    end

  end


% -----------------------------------------------------------------------
function pass = ping_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd)

  nnq = size(qd.qPhi, 2);
  nnv = size(qd.vPhi, 2);
  nnu = size(qd.uPhi, 2);
  nnw = size(qd.wPhi, 2);

  ep = 1e-1;
  tol = 1e-12;

  [Rq0, Ru0, Rl0, Rq_Q0, Rq_U0, Rq_L0, Ru_Q0, Ru_U0, Ru_L0, Rl_Q0, Rl_U0, Rl_L0] = ...
  hdg_residual_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

  pass = true;

  % loop over q dof
  for j=1:nnq
    QE(j) = QE(j) + ep;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

    QE(j) = QE(j) - ep;

    diff = [Rq-Rq0; Ru-Ru0; Rl-Rl0]/ep;
    exact = 0.5*[Rq_Q(:,j)+Rq_Q0(:,j); Ru_Q(:,j)+Ru_Q0(:,j); Rl_Q(:,j)+Rl_Q0(:,j)];
    err = abs(diff - exact);
    if any(err > tol)
      err
      warning('Ping error: Q residuals, j=%d', j)
      pass = false;
      return;
    end
  end

  % loop over u dof
  for j=1:nnu
    UE(j) = UE(j) + ep;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

    UE(j) = UE(j) - ep;

    diff = [Rq-Rq0; Ru-Ru0; Rl-Rl0]/ep;
    exact = 0.5*[Rq_U(:,j)+Rq_U0(:,j); Ru_U(:,j)+Ru_U0(:,j); Rl_U(:,j)+Rl_U0(:,j)];
    err = abs(diff - exact);
    if any(err > tol)
      err
      warning('Ping error: U residuals, j=%d', j)
      pass = false;
      return;
    end
  end

  % loop over l dof
  for j=1:2
    LE(j) = LE(j) + ep;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    hdg_residual_elem(QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

    LE(j) = LE(j) - ep;

    diff = [Rq-Rq0; Ru-Ru0; Rl-Rl0]/ep;
    exact = 0.5*[Rq_L(:,j)+Rq_L0(:,j); Ru_L(:,j)+Ru_L0(:,j); Rl_L(:,j)+Rl_L0(:,j)];
    err = abs(diff - exact);
    if any(err > tol)
      err
      warning('Ping error: L residual j=%d', j)
      pass = false;
      return;
    end
  end


