function [dQ, dU, dL] = hdgdpg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU, dL] = hdgdpg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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
  %   qd : quadrature data with optimal test functions [struct]
  %
  % OUTPUTS:
  %   d{Q,U} : update for {Q,U} [nelem*nn{q,u}]
  %   dL : update for L [nelem-1]

  % 1. Build system (K*L=F) (F=-Residual)
  [K, F] = system(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

  % 2. Solve reduced HDG system (linear solve)
  dL = K\F;

  % 3. Post-linear solve
  [dQ, dU] = post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd);


% -----------------------------------------------------------------------
function [K, R] = system(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [K, R] = system(Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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
  
  xg = qd.x*md.dx;

  [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
  residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				 Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  if md.ne ~= 1
    B(1,2) = BK(2,2);
    R(1) = BR(2);
  end

  lisb = false;
  risb = false;

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);

    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

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
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				   Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    B(md.ne-1,2) = B(md.ne-1,2) + BK(1,1);
    R(md.ne-1) = R(md.ne-1) + BR(1);
  end

  K = spdiags(B, -1:1, md.ne-1, md.ne-1);


% -----------------------------------------------------------------------
function [dQ, dU] = post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
  % function [dQ, dU] = post(dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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
  residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

  [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			    Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

  dQ(1:nnq) = dQE;
  dU(1:nnu) = dUE;


  lisb = false;
  risb = false;

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);
    dLE = dL(elem-1:elem);

    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

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
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			      Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((md.ne-1)*nnq+1:md.ne*nnq) = dQE;
    dU((md.ne-1)*nnu+1:md.ne*nnu) = dUE;
  end


% -----------------------------------------------------------------------
function [BK, BR] = static_condensation(Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
					Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % initialize matrix and vector
  % BK = D
  % BR = H
  BK = Rl_L;
  BR = -Rl;

  A = build_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  if q_present
    BK = BK - [Rl_Q, Rl_U]*(A\(Rq_L+Ru_L));
  else
    BK = BK - Rl_U*(A\Ru_L);
  end

  if q_present
    BR = BR - [Rl_Q, Rl_U]*(A\(-(Rq+Ru)));
  else
    BR = BR - Rl_U*(A\(-Ru));
  end


% -----------------------------------------------------------------------
function [dQ, dU] = qu_backsolve(dL, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				 Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  nnq = size(Ru_Q, 2);
  nnu = size(Ru_U, 2);

  A = build_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  % F - B*dL
  F = -(Rq+Ru);
  dQU = A\(F - (Rq_L+Ru_L)*dL);

  if q_present
    dQ = dQU(1:nnq);
    dU = dQU(nnq+(1:nnu));
  else
    dQ = zeros([nnq,1]);
    dU = dQU(1:nnu);
  end


% -----------------------------------------------------------------------
function A = build_A(Rq_Q, Rq_U, Ru_Q, Ru_U, q_present)

  nn = size(Ru_U, 1);
  nnq = size(Ru_Q, 2);
  nnu = size(Ru_U, 2);

  A = zeros([nn, nn]);

  if q_present
    A(:,1:nnq) = Rq_Q;
    A(:,nnq+(1:nnu)) = Rq_U;
  else
    nnq = 0;
  end
  A(:,1:nnq) = A(:,1:nnq) + Ru_Q;
  A(:,nnq+(1:nnu)) = A(:,nnq+(1:nnu)) + Ru_U;
