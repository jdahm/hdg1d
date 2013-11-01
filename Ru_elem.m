function [R, R_Q, R_U, R_UH] = Ru_elem(Q, U, uh, xg, td, fd, sd, dx, lisb, risb, qd)
  % function [R, R_Q, R_U, R_UH] = Ru_elem(Q, U, uh, xg, td, fd, sd, dx, qd)
  %
  % PURPOSE: Computes the bilinear form R((q,u,uh),w) : Qh x Uh x Mh x Wh -> Wh*
  % restricted to an element and linearization wrt each input. This is:
  %
  % (u_t,w) - (H,grad(w)) + <HH,w> + (S,w) = 0 forall w in Wh
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nn{q,u}]
  %   uh : u|{0,1}, that is, the trace of u at the boundaries of the element
  %   xg : global positions of quadrature nodes [nq]
  %   td : time data [struct]
  %   dx : size of the element [1]
  %   qd : quadrature data [struct]
  %   fd : flux data [struct]
  %   sd : source data [struct]
  %
  % OUTPUTS:
  %   R : bilinear form [nnw]
  %   R_{Q,U,UH} : linearization of R wrt {q,u,uh}
  %

  nnq = size(qd.qPhi, 2);
  nnu = size(qd.uPhi, 2);

  q = qd.qPhi * Q;
  u = qd.uPhi * U;

  q0 = qd.qPhi0 * Q;
  q1 = qd.qPhi1 * Q;
  u0 = qd.uPhi0 * U;
  u1 = qd.uPhi1 * U;

  n0 = -1.;
  n1 = 1.;

  [h0, h_q0, h_UH0] = flux(q0, uh(1), fd);
  h_u0 = 0.;
  % dot with normal
  h0 = h0*n0;
  h_q0 = h_q0*n0;
  h_UH0 = h_UH0*n0;

  [s0, s_u0, s_UH0] = flux_stab(u0, uh(1), n0, fd);

  f0 = h0 + s0;
  f_Q0 = h_q0*qd.qPhi0;
  f_U0 = (h_u0+s_u0)*qd.uPhi0;
  f_UH0 = h_UH0 + s_UH0;

  [h1, h_q1, h_UH1] = flux(q1, uh(2), fd);
  h_u1 = 0.;
  % dot with normal
  h1 = h1*n1;
  h_q1 = h_q1*n1;
  h_UH1 = h_UH1*n1;

  [s1, s_u1, s_UH1] = flux_stab(u1, uh(2), n1, fd);

  f1 = h1 + s1;
  f_Q1 = h_q1*qd.qPhi1;
  f_U1 = (h_u1+s_u1)*qd.uPhi1;
  f_UH1 = h_UH1 + s_UH1;


  [f, f_q, f_u] = flux(q, u, fd);
  f_Q = repmat(f_q, [1, nnq]) .* qd.qPhi;
  f_U = repmat(f_u, [1, nnu]) .* qd.uPhi;

  [s, s_q, s_u] = source(q, u, xg, fd, sd);
  s_Q = repmat(s_q, [1, nnq]) .* qd.qPhi;
  s_U = repmat(s_u, [1, nnu]) .* qd.uPhi;

  wgPhi = qd.wGPhi/dx;

  % Ru = (u_t,w) - (H,grad(w)) + <HH,w> + (S,w)

  nq = size(qd.wPhi,1);

  % I1 = (u_t,w)
  I1   = dx * qd.wPhi'*qd.dw*ones([nq,1])*td.u_t;
  I1_U = dx * qd.wPhi'*qd.dw*(td.c*qd.uPhi);

  % I2 = (H,grad(w)) = (H,w_x) (1D)
  I2   = dx * wgPhi'*qd.dw*f;
  I2_Q = dx * wgPhi'*qd.dw*f_Q;
  I2_U = dx * wgPhi'*qd.dw*f_U;

  % I3 = <HH,w> = f1*wPhi1' + f0*wPhi0' (1D)
  I3 = f1*qd.wPhi1' + f0*qd.wPhi0';
  I3_Q = qd.wPhi1'*f_Q1 + qd.wPhi0'*f_Q0;
  I3_U = qd.wPhi1'*f_U1 + qd.wPhi0'*f_U0;
  I3_UH(:,1) = qd.wPhi0'*f_UH0;
  I3_UH(:,2) = qd.wPhi1'*f_UH1;

  % I4 = (S,w)
  I4   = dx * qd.wPhi'*qd.dw*s;
  I4_Q = dx * qd.wPhi'*qd.dw*s_Q;
  I4_U = dx * qd.wPhi'*qd.dw*s_U;

  R = I1 - I2 + I3 + I4;
  R_Q = -I2_Q + I3_Q + I4_Q;
  R_U = I1_U - I2_U + I3_U + I4_U;
  R_UH = I3_UH;


%{
%This is the version with the adjoint consistency fix, but since we
%aren't using a Riemann solver on the boundaries, this code is
%deprecated
  if lisb
    [h0, h_q0, h_u0, h_UH0] = flux_bc(q0, u0, uh(1), n0, fd);
    if strcmp(scheme, 'dpg')
      % adjoint consistency fix
      taut = 0.5*(fd.a*n0+abs(fd.a*n0)) + fd.b/fd.vl;
      c0 = -0.5*(fd.a*n0+abs(fd.a*n0))/taut;
    else
      c0 = 1.0;
    end
  else
    [h0, h_q0, h_UH0] = flux(q0, uh(1), fd);
    h_u0 = 0.;
    % dot with normal
    h0 = h0*n0;
    h_q0 = h_q0*n0;
    h_UH0 = h_UH0*n0;
    if strcmp(scheme, 'dpg')
      c0 = fd.c;
    else
      c0 = 1.0;
    end
  end
  [s0, s_u0, s_UH0] = flux_stab(u0, uh(1), n0, fd);
  f0 = h0 + c0*s0;
  f_Q0 = h_q0*qd.qPhi0;
  f_U0 = (h_u0+c0*s_u0)*qd.uPhi0;
  f_UH0 = h_UH0 + c0*s_UH0;

  if risb
    [h1, h_q1, h_u1, h_UH1] = flux_bc(q1, u1, uh(2), n1, fd);
    if strcmp(scheme, 'dpg')
      taut = 0.5*(fd.a*n1+abs(fd.a*n1)) + fd.b/fd.vl;
      c1 = -0.5*(fd.a*n1+abs(fd.a*n1))/taut;
    else
      c1 = 1.0;
    end
  else
    [h1, h_q1, h_UH1] = flux(q1, uh(2), fd);
    h_u1 = 0.;
    % dot with normal
    h1 = h1*n1;
    h_q1 = h_q1*n1;
    h_UH1 = h_UH1*n1;
    if strcmp(scheme, 'dpg')
      c1 = fd.c;
    else
      c1 = 1.0;
    end
  end
  [s1, s_u1, s_UH1] = flux_stab(u1, uh(2), n1, fd);
  f1 = h1 + c1*s1;
  f_Q1 = h_q1*qd.qPhi1;
  f_U1 = (h_u1+c1*s_u1)*qd.uPhi1;
  f_UH1 = h_UH1 + c1*s_UH1;
%}
