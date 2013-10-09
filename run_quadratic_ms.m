% Quadratic manufactured solution

fd.a = 1;
fd.b = 1;
fd.q_present = true;
fd.stab_type = 'centered';
fd.vl = 1.0;

sd.present = true;
sd.type = 'ms_quadratic';

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 7;
pv = 7;
pu = 7;
pw = 7;

md = mesh(0., 2., 2);

xnq = create_nodes(pq, 'SegLagrangeGauss');
xnv = create_nodes(pv, 'SegLagrangeGauss');
xnu = create_nodes(pu, 'SegLagrangeGauss');
xnw = create_nodes(pw, 'SegLagrangeGauss');

lbd.type = 'd';
lbd.data = 0.;
rbd.type = 'd';
rbd.data = 0.;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1);
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);
