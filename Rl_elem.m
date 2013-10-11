function [R, R_Q, R_U, R_UH] = Rl_elem(Q, U, uh, fd, lisb, risb, qd)
  % function [R, R_Q, R_U, R_UH] = Rl_elem(Q, U, uh0, uh1, td, fd, sd, dx, qd)
  %
  % PURPOSE: Computes the bilinear form R((q,u,uh),m) : Qh x Uh x Mh x Mh -> Mh*
  % restricted to an element and linearization wrt each input. This is:
  %
  % <[HH],m> = 0 forall m in Mh
  %
  % NOTE: On interior elements, computes a one-sided stabilized flux on both edges,
  % on boundary edges, computes the flux H(q,uh)
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nn{q,u}]
  %   uh : u|{0,1}, that is, the trace of u at the boundaries of the element
  %   fd : flux data [struct]
  %   lisb : true if left is boundary, else false [bool]
  %   risb : true if right is boundary, else false [bool]
  %   qd : quadrature data [struct]
  %
  % OUTPUTS:
  %   R : bilinear form [1]
  %   R_{Q,U,UH} : linearization of R wrt {q,u,uh}
  %

  nnq = size(qd.qPhi, 2);
  nnu = size(qd.uPhi, 2);

  q0 = qd.qPhi0 * Q;
  q1 = qd.qPhi1 * Q;
  u0 = qd.uPhi0 * U;
  u1 = qd.uPhi1 * U;

  n0 = -1.;
  n1 = 1.;

  if lisb
    [h0, h_q0, h_u0, h_UH0] = flux_bc(q0, u0, uh(1), n0, fd);
  else
    [h0, h_q0, h_UH0] = flux(q0, uh(1), fd);
    h_u0 = 0.;
    % dot with normal
    h0 = h0*n0;
    h_q0 = h_q0*n0;
    h_UH0 = h_UH0*n0;
  end
  [s0, s_u0, s_UH0] = flux_stab(u0, uh(1), n0, fd);
  f0 = h0 + s0;
  f_Q0 = h_q0*qd.qPhi0;
  f_U0 = (h_u0+s_u0)*qd.uPhi0;
  f_UH0 = h_UH0 + s_UH0;

  if risb
    [h1, h_q1, h_u1, h_UH1] = flux_bc(q1, u1, uh(2), n1, fd);
  else
    [h1, h_q1, h_UH1] = flux(q1, uh(2), fd);
    h_u1 = 0.;
    % dot with normal
    h1 = h1*n1;
    h_q1 = h_q1*n1;
    h_UH1 = h_UH1*n1;
  end
  [s1, s_u1, s_UH1] = flux_stab(u1, uh(2), n1, fd);
  f1 = h1 + s1;
  f_Q1 = h_q1*qd.qPhi1;
  f_U1 = (h_u1+s_u1)*qd.uPhi1;
  f_UH1 = h_UH1 + s_UH1;

  % Rl = <[HH],m>
  R(1,1)    = f0;
  R(2,1)    = f1;
  R_Q(1,:)  = f_Q0;
  R_Q(2,:)  = f_Q1;
  R_U(1,:)  = f_U0;
  R_U(2,:)  = f_U1;
  R_UH(1,1) = f_UH0;
  R_UH(2,2) = f_UH1;
