function [R, R_Q, R_U, R_UH] = Ru_elem(Q, U, uh, td, fd, sd, dx, qd)
  % function [R, R_Q, R_U, R_UH] = Ru_elem(Q, U, uh, td, fd, sd, dx, qd)
  %
  % PURPOSE: Computes the bilinear form Ru((q,u,uh),w) : Qh x Uh x Mh x Wh -> Wh*
  % restricted to an element and linearization wrt each input. This is:
  %
  % (u_t,w) - (H,grad(w)) + <HH,w> + (S,w) = 0 forall w in Wh
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nn{q,u}]
  %   uh : u|{0,1}, that is, the trace of u at the boundaries of the element
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

  q = qd.qPhi * Q;
  u = qd.uPhi * U;

  q0 = qd.qPhi0 * Q;
  q1 = qd.qPhi1 * Q;
  u0 = qd.uPhi0 * U;
  u1 = qd.uPhi1 * U;

  [f0, f_q0, f_u0, f_UH0] = flux_oneside(q0, u0, uh(1), fd);
  f_Q0 = bsxfun(@times, f_q0, qd.qPhi0);
  f_U0 = bsxfun(@times, f_u0, qd.uPhi0);

  [f1, f_q1, f_u1, f_UH1] = flux_oneside(q1, u1, uh(2), fd);
  f_Q1 = bsxfun(@times, f_q1, qd.qPhi1);
  f_U1 = bsxfun(@times, f_u1, qd.uPhi1);

  [f, f_q, f_u] = flux(q,u,fd);
  f_Q = bsxfun(@times, f_q, qd.qPhi);
  f_U = bsxfun(@times, f_u, qd.uPhi);

  [s, s_q, s_u] = source(q,u,sd);
  s_Q = bsxfun(@times, s_q, qd.qPhi);
  s_U = bsxfun(@times, s_u, qd.uPhi);

  wgPhi = qd.wGPhi*dx;

  % Ru = (u_t,w) - (H,grad(w)) + <HH,w> + (S,w)

  nq = size(qd.wPhi,1);

  % I1 = (u_t,w)
  I1   = dx * qd.wPhi'*qd.dw*ones([nq,1])*td.u_t;
  I1_U = dx * qd.wPhi'*qd.dw*(td.c*qd.uPhi);

  % I2 = (H,grad(w)) = (H,w_x) (1D)
  I2   = dx * wgPhi'*qd.dw*f;
  I2_Q = dx * wgPhi'*qd.dw*f_Q;
  I2_U = dx * wgPhi'*qd.dw*f_U;

  % I3 = <HH,w> = f1*wPhi1' - f0*wPhi0' (1D)
  I3 = f1*qd.wPhi1' - f0*qd.wPhi0';
  I3_Q = bsxfun(@times, qd.wPhi1', f_Q1) - bsxfun(@times, qd.wPhi0', f_Q0);
  I3_U = bsxfun(@times, qd.wPhi1', f_U1) - bsxfun(@times, qd.wPhi0', f_U0);
  I3_UH(:,1) = -bsxfun(@times, qd.wPhi0', f_UH0);
  I3_UH(:,2) = bsxfun(@times, qd.wPhi1', f_UH1);

  % I4 = (S,w)
  I4   = dx * qd.wPhi'*qd.dw*s;
  I4_Q = dx * qd.wPhi'*qd.dw*s_Q;
  I4_U = dx * qd.wPhi'*qd.dw*s_U;

  R = I1 - I2 + I3 + I4;
  R_Q = -I2_Q + I3_Q + I4_Q;
  R_U = I1_U - I2_U + I3_U + I4_U;
  R_UH = I3_UH;
