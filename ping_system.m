function pass = ping_system(scheme, Q, U, L, lbd, rbd, md, td, fd, sd, qd)

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

  xg = qd.x*md.dx;

  pass_elem = ping_elem(scheme, QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
  if ~pass_elem
     pass = false;
     return;
  end

  lisb = false;
  risb = false;

  % interior elements
  for elem=2:md.ne-1
    QE = Q((elem-1)*nnq+1:elem*nnq);
    UE = U((elem-1)*nnu+1:elem*nnu);
    LE = L(elem-1:elem);

    xg = (md.ne-1)*md.dx+qd.x*md.dx;

    pass_elem = ping_elem(scheme, QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
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

    pass_elem = ping_elem(scheme, QE, UE, LE, xg, lisb, risb, lbd, rbd, md.dx, td, fd, sd, qd);
    if ~pass_elem
      pass = false;
      return;
    end

  end


% -----------------------------------------------------------------------
function pass = ping_elem(scheme, QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd)

  nnq = size(qd.qPhi, 2);
  nnv = size(qd.vPhi, 2);
  nnu = size(qd.uPhi, 2);
  nnw = size(qd.wPhi, 2);

  ep = 1e-1;
  tol = 1e-12;

  [Rq0, Ru0, Rl0, Rq_Q0, Rq_U0, Rq_L0, Ru_Q0, Ru_U0, Ru_L0, Rl_Q0, Rl_U0, Rl_L0] = ...
  residual_elem(scheme, QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

  pass = true;

  % loop over q dof
  for j=1:nnq
    QE(j) = QE(j) + ep;

    [Rq, Ru, Rl, Rq_Q, Rq_U, Rq_L, Ru_Q, Ru_U, Ru_L, Rl_Q, Rl_U, Rl_L] = ...
    residual_elem(scheme, QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

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
    residual_elem(scheme, QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

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
    residual_elem(scheme, QE, UE, LE, xglob, lisb, risb, lbd, rbd, dx, td, fd, sd, qd);

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
