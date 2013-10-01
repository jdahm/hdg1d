function [R, R_Q, R_U, R_UH] = Rq_elem(Q, U, uh, dx, qd)
  % function [R, R_Q, R_U, R_UH] = Rq_elem(Q, U, uh, dx, qd)
  %
  % PURPOSE: Computes the bilinear form Rq((q,u,uh),v) : Qh x Uh x Mh x Vh -> Vh*
  % restricted to an element and linearization wrt each input. This is:
  %
  % (q,v) + (u,div(v)) - <uh,v.n> = 0 forall v in Vh
  %
  % INPUTS:
  %   {Q,U} : basis coefficients for {grad(u),u} [nn{q,u}]
  %   uh : u|{0,1}, that is, the trace of u at the boundaries of the element
  %   dx : size of the element [1]
  %   qd : quadrature data [struct]
  %
  % OUTPUTS:
  %   R : bilinear form [nnv]
  %   R_{Q,U,UH} : linearization of R wrt {q,u,uh}
  %

  q = qd.qPhi * Q;
  u = qd.uPhi * U;

  vgPhi = qd.vGPhi*dx;

  % Rq = (q,v) + (u,div(v)) - <uh,v.n>

  % I1 = (q,v)
  I1   = dx * qd.vPhi'*qd.dw*q;
  I1_Q = dx * qd.vPhi'*qd.dw*qd.qPhi;

  % I2 = (u,div(v)) = (u,v_x) (1D)
  I2   = dx * vgPhi'*qd.dw*u;
  I2_U = dx * vgPhi'*qd.dw*qd.uPhi;

  % I3 = <uh,v.n> = uh1*vPhi1' - uh0*vPhi0' (1D)
  I3 = uh(2)*qd.vPhi1' - uh(1)*qd.vPhi0';
  I3_UH(:,1) = -1.0*qd.vPhi0';
  I3_UH(:,2) = qd.vPhi1';

  R = I1 + I2 - I3;
  R_Q = I1_Q;
  R_U = I2_U;
  R_UH = -I3_UH;
