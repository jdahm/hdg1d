% Example Advection-diffusion HDG solution

fd.a = 1.0;
fd.b = 1.0;
fd.q_present = true;
fd.stab = 1.0;

sd.a = 1.0;
sd.b = 0.0;

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 3;
pv = 3;
pu = 3;
pw = 3;

md = mesh(0.0, 1.0, 4);

xnq = chebyshev_roots(pq+1)';
xnv = chebyshev_roots(pv+1)';
xnu = chebyshev_roots(pu+1)';
xnw = chebyshev_roots(pw+1)';

lbd.type = 'd';
lbd.data = 1.0;
rbd.type = 'd';
rbd.data = 1.0;

[Q, U, L] = initialize(pq, pu, ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

figure(1);
h1 = plot_elems(md.xs, md.xe, xnu, U, 10);
%plot_trace(xs, xe, b0d, b1d, L);
%figure(2);
%h2 = plot_elems(xs, xe, xnq, Q, 10);
