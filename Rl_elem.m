function [R, R_Q, R_U, R_UH] = Rl_elem(Q, U, uh, fd, qd)
  % function [R, R_Q, R_U, R_UH] = Ru_elem(Q, U, uh0, uh1, td, fd, sd, dx, qd)
  %
  % PURPOSE: Computes the bilinear form Rl((q,u,uh),w) : Qh x Uh x Mh x Mh -> Mh*
  % restricted to an element and linearization wrt each input. This is:
  %
  % <[HH],m> = 0 forall m in Mh
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nn{q,u}]
  %   uh : u|{0,1}, that is, the trace of u at the boundaries of the element
  %   qd : quadrature data [struct]
  %
  % OUTPUTS:
  %   R : bilinear form [1]
  %   R_{Q,U,UH} : linearization of R wrt {q,u,uh}
  %

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

  % Rl = <[HH],m>

  R(1,1)    = f0;
  R(2,1)    = f1;
  R_Q(1,:)  = f_Q0;
  R_Q(2,:)  = f_Q1;
  R_U(1,:)  = f_U0;
  R_U(2,:)  = f_U1;
  R_UH(1,1) = f_UH0;
  R_UH(2,2) = f_UH1;
