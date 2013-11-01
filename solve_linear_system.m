function [dQ, dU, dL] = solve_linear_system(scheme, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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
  [K, F] = system(scheme, Q, U, L, lbd, rbd, md, td, fd, sd, qd);

  % 2. Solve reduced HDG system (linear solve)
  dL = K\F;

  % 3. Post-linear solve
  [dQ, dU] = post(scheme, dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd);


% -----------------------------------------------------------------------
function [K, R] = system(scheme, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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

  % loop over elements
  for elem=1:md.ne
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = [0.0; 0.0];
    lisb = true;
    risb = true;

    if elem ~= 1
       LE(1) = L(elem-1);
       lisb = false;
    end
    if elem ~= md.ne
       LE(2) = L(elem);
       risb = false;
    end

    % global positions of quad points
    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [BK, BR] = static_condensation(scheme, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				   Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    if elem ~= 1
      % self-block
      B(elem-1,2) = B(elem-1,2) + BK(1,1);
      % off-diagonal block
      B(elem-1,1) = B(elem-1,1) + BK(2,1);
      R(elem-1) = R(elem-1) + BR(1);
    end
    if elem ~= md.ne
      B(elem,2)   = B(elem,2)   + BK(2,2);
      % off-diagonal entries
      B(elem,3) = B(elem,3) + BK(1,2);
      R(elem) = R(elem) + BR(2);
    end
  end

  K = spdiags(B, -1:1, md.ne-1, md.ne-1);


% -----------------------------------------------------------------------
function [dQ, dU] = post(scheme, dL, Q, U, L, lbd, rbd, md, td, fd, sd, qd)
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

  % loop over elements
  for elem=1:md.ne
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = [0.0; 0.0];
    dLE = [0.0; 0.0];
    lisb = true;
    risb = true;

    if elem ~= 1
       LE(1) = L(elem-1);
       dLE(1) = dL(elem-1);
       lisb = false;
    end
    if elem ~= md.ne
       LE(2) = L(elem);
       dLE(2) = dL(elem);
       risb = false;
    end

    % global positions of quad points
    xg = (elem-1)*md.dx+qd.x*md.dx;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    residual_elem(QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);

    [dQE, dUE] = qu_backsolve(scheme, dLE, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
			      Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, fd.q_present);

    dQ((elem-1)*nnq+1:elem*nnq) = dQE;
    dU((elem-1)*nnu+1:elem*nnu) = dUE;
  end


% -----------------------------------------------------------------------
function [BK, BR] = static_condensation(scheme, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
					Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % initialize matrix and vector
  % BK = D
  % BR = H
  BK = Rl_L;
  BR = -Rl;

  A = build_and_preprocess_A(scheme, Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  if strcmp(scheme, 'hdg')

    % subtract C(A^{-1}B) from BK
    [iABQ, iABU] = hdg_apply_Ainv(A, Rq_L, Ru_L, q_present);
    BK = BK - Rl_U*iABU;
    if q_present
      BK = BK - Rl_Q*iABQ;
    end

    % subtract C(A^{-1}F) from BR, where F=-[Rq;Ru]
    [iAFQ, iAFU] = hdg_apply_Ainv(A, -Rq, -Ru, q_present);
    BR = BR - Rl_U*iAFU;
    if q_present
      BR = BR - Rl_Q*iAFQ;
    end

  elseif strcmp(scheme, 'dpg')

    % subtract C(A^{-1}B) from BK
    if q_present
      BK = BK - [Rl_Q, Rl_U]*(A\(Rq_L+Ru_L));
    else
      BK = BK - Rl_U*(A\Ru_L);
    end

    % subtract C(A^{-1}F) from BR, where F=-(Rq+Ru)
    if q_present
      BR = BR - [Rl_Q, Rl_U]*(A\(-(Rq+Ru)));
    else
      BR = BR - Rl_U*(A\(-Ru));
    end

  end


% -----------------------------------------------------------------------
function [dQ, dU] = qu_backsolve(scheme, dL, Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, ...
				 Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L, q_present)

  % pre-process A matrix
  A = build_and_preprocess_A(scheme, Rq_Q, Rq_U, Ru_Q, Ru_U, q_present);

  if strcmp(scheme, 'hdg')

    % R = F - B*dL
    Rq = -Rq - Rq_L*dL;
    Ru = -Ru - Ru_L*dL;

    % QU = A^{-1}*R
    [dQ, dU] = hdg_apply_Ainv(A, Rq, Ru, q_present);

  elseif strcmp(scheme, 'dpg')

    nnq = size(Ru_Q, 2);
    nnu = size(Ru_U, 2);

    % F - B*dL
    F = -(Rq+Ru);

    % QU = A^{-1}*R
    dQU = A\(F - (Rq_L+Ru_L)*dL);

    if q_present
      dQ = dQU(1:nnq);
      dU = dQU(nnq+(1:nnu));
    else
      dQ = zeros([nnq,1]);
      dU = dQU(1:nnu);
    end

  end


% -----------------------------------------------------------------------
function A = build_and_preprocess_A(scheme, Rq_Q, Rq_U, Ru_Q, Ru_U, q_present)
  % function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)
  %
  % PURPOSE: Takes a block matrix A
  %   A = [AQQ, AQU; AUQ, AUU]  if q_present is true
  %   A = [AUU]                 otherwise
  % and prepares it to be quickly inverted.
  %
  % INPUTS:
  %   scheme : scheme name [string]
  %   {AQQ, AQU, AUQ, AUU} : matrix blocks
  %   q_present : indicator of number of blocks
  %

  if strcmp(scheme, 'hdg')

    % OUTPUTS:
    %   A : [struct]
    %   A.QQ : stays same (A.QQ)
    %   A.QU : becomes A.QQ^{-1}*A.UQ
    %   A.UQ : K^{-1}*A.UQ where K=(A.UU-A.UQ*A.QQ^{-1}*A.UQ)
    %   A.UU : K^{-1}

    A.QQ = Rq_Q;
    A.QU = Rq_U;
    A.UQ = Ru_Q;
    A.UU = Ru_U;

    if q_present % needed because AQQ should not be invertible unless q is present
      % AQU = AQQ^{-1}*AQU
      A.QU = A.QQ\A.QU;
      % AUU -= AUQ*AQU
      A.UU = A.UU-A.UQ*A.QU;
      % AUQ = AUU^{-1}*AUQ
      A.UQ = A.UU\A.UQ;
    end

  elseif strcmp(scheme, 'dpg')

    % OUTPUTS:
    %   A : matrix without blocks

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

  else

    error('unknown method');

  end


% -----------------------------------------------------------------------
function [Q, U] = hdg_apply_Ainv(A, Q, U, q_present)
  % function [Q, U] = hdg_apply_Ainv(A, Q, U, q_present)
  %
  % PURPOSE: Takes a pre-processed block matrix A and performs
  %   A^{-1}*[Q;U]
  % where Q and U can either form a matrix or vector.
  %
  % INPUTS:
  %   A.{AQQ, AQU, AUQ, AUU} : matrix blocks after pre-processing
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
    Q = A.QQ\Q;
  end
  % W = -AUU^{-1}*U
  W = -A.UU\U;
  if q_present
    % W += AUQ*Q
    W = W+A.UQ*Q;
    % Q += AQU*W
    Q = Q+A.QU*W;
  else
    Q = zeros(size(Q));
  end
  % U = -W
  U = -W;
