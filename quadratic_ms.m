% Quadratic manufactured solution

fd.a = 1.0;
fd.b = 1.0;
fd.q_present = true;
fd.stab = 1.0;

sd.present = true;
sd.type = 'ms_quadratic';

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 2;
pv = 2;
pu = 2;
pw = 2;

md = mesh(0., 2., 5);

xnq = linspace(0,1,pq+1)';
xnv = linspace(0,1,pv+1)';
xnu = linspace(0,1,pu+1)';
xnw = linspace(0,1,pw+1)';

lbd.type = 'd';
lbd.data = 0.0;
rbd.type = 'd';
rbd.data = 0.0;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd)

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1);
h1 = plot_elems(md.xs, md.xe, xnu, U, 10);
%plot_trace(xs, xe, b0d, b1d, L);
%figure(2);
%h2 = plot_elems(xs, xe, xnq, Q, 10);
