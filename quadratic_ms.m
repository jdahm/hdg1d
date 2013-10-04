% Quadratic manufactured solution

fd.a = 1;
fd.b = 1e-3;
fd.q_present = false;
fd.vl = 1e-3*fd.b;

sd.present = false;
sd.type = 'ms_quadratic';

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 5;
pv = 5;
pu = 5;
pw = 5;

md = mesh(0., 2., 10);

xnq = linspace(0,1,pq+1)';
xnv = linspace(0,1,pv+1)';
xnu = linspace(0,1,pu+1)';
xnw = linspace(0,1,pw+1)';

lbd.type = 'd';
lbd.data = 2;
rbd.type = 'n';
rbd.data = 0;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1);
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);
%plot_trace(xs, xe, b0d, b1d, L);
%figure(2);
%h2 = plot_elems(xs, xe, xnq, Q, 10);
