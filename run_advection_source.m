% Advection with a source term

scheme = 'hdg';
%scheme = 'dpg';

fd.a = 1;
fd.b = 0;
fd.q_present = false;
fd.stab_type = 'upwind';
fd.vl = 1.0;
fd.c = 1.0;

sd.present = true;
sd.type = 'linear';
sd.a = -5;
sd.b = 0;

td.c = 0.0;
td.u_t = 0.0;

pq = 2;
pu = 2;

pv = 2;
pw = 2;

qbasis = 'SegLagrange';
ubasis = 'SegLagrange';

% optimal test function weights (higher weight means more emphasis on
% boundary portion of output)
testd.wleft = 1e3;
testd.wright = 1e3;

md = mesh(0., 1., 32);

lbd.type = 'd';
lbd.data = 1.;
rbd.type = 'n';
%rbd.data = 0.;

xnq = create_nodes(pq, qbasis);
xnv = create_nodes(pv, qbasis);
xnu = create_nodes(pu, ubasis);
xnw = create_nodes(pw, ubasis);

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quad_data(xnq, xnu, xnv, xnw);

if strcmp(scheme, 'dpg')

  % compute optimal test functions
  Vopt = compute_optimal_test_functions(fd, qd, md, sd, lbd, rbd, testd);

  %% plot test functions
  figure(3); clf;
  plot_test_functions(Vopt.w, xnw, qd, 1000);

  % fill DPG quadrature data
  qd = dpg_quad_data(qd, Vopt.v, Vopt.w);
end

[dQ, dU, dL] = solve_linear_system(scheme, Q, U, L, lbd, rbd, ...
				   md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 1000);
if md.ne > 1
  h1 = plot_traces(md.xs, md.xe, L);
end

% plot exact solution
x = linspace(0,1,1000);
y = exp(-sd.a/fd.a*x);
hold on;
plot(x,y,'--');
hold off;

erru = error_norm(Q, U, L, xnq, xnu, pu+1, md, fd, ...
		  'L2', @(fd, x) exp(5/fd.a*x))
errup = error_norm(Q, PU, L, xnq, pxnu, pu+2, md, fd, ...
		   'L2', @(fd, x) exp(5/fd.a*x))
