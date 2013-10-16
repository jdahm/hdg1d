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
function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)
  % function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)
  %
  % PURPOSE: Takes a block matrix A
  %   A = [AQQ, AQU; AUQ, AUU]  if q_present is true
  %   A = [AUU]                 otherwise
  % and prepares it to be quickly inverted.
  %
  % INPUTS:
  %   {AQQ, AQU, AUQ, AUU} : matrix blocks
  %   q_present : indicator of number of blocks
  %
  % OUTPUTS:
  %   AQQ : stays same (AQQ)
  %   AQU : becomes AQQ^{-1}*AUQ
  %   AUQ : K^{-1}*AUQ where K=(AUU-AUQ*AQQ^{-1}*AUQ)
  %   AUU : K^{-1}
  %

  if q_present % needed because AQQ should not be invertible unless q is present
    % AQU = AQQ^{-1}*AQU
    AQU = AQQ\AQU;
    % AUU -= AUQ*AQU
    AUU = AUU-AUQ*AQU;
    % AUQ = AUU^{-1}*AUQ
    AUQ = AUU\AUQ;
  end


% -----------------------------------------------------------------------
function [Q, U] = apply_Ainv(AQQ, AQU, AUQ, AUU, Q, U, q_present)
  % function [Q, U] = apply_Ainv(AQQ, AQU, AUQ, AUU, Q, U, q_present)
  %
  % PURPOSE: Takes a pre-processed block matrix A and performs
  %   A^{-1}*[Q;U]
  % where Q and U can either form a matrix or vector.
  %
  % INPUTS:
  %   {AQQ, AQU, AUQ, AUU} : matrix blocks after pre-processing
  %   Q : block [nq, c]
  %   U : block [nu, c]
  %   q_present : indicator of number of blocks
  %
  % OUTPUTS:
  %   Q : first nq rows of A^{-1}*[Q;U]
  %   U : last nu rows of A^{-1}*[Q;U]
  %

  if q_present
    % needed since AQQ should not be invertible unless q is present
    % Q = AQQ^{-1}*Q
    Q = AQQ\Q;
  end
  % W = -AUU^{-1}*U
  W = -AUU\U;
  if q_present
    % W += AUQ*Q
    W = W+AUQ*Q;
    % Q += AQU*W
    Q = Q+AQU*W;
  else
    Q = zeros(size(Q));
  end
  % U = -W
  U = -W;

