% Boundary layer
clear all; close all; format longe; clc;

scheme = 'hdg';
scheme = 'dpg';

fd.a = 1;
fd.b = 1e-1;
fd.q_present = true;
fd.stab_type = 'centered';
fd.vl = 1.0;
fd.c = 1.0;
fd.c = 0.0;

sd.present = true;
sd.type = 'linear';
sd.a = -3;
sd.b = 0;

td.c = 0.0;
td.u_t = 0.0;

pq = 1;
pu = 1;

pv = 20;
pw = 20;

% optimal test function weights (higher weight means more emphasis on
% boundary portion of output)
testd.wleft = 1e12;
testd.wright = 1e12;

md = mesh(0., 1., 5);

qbasis = 'SegLagrangeGauss';
ubasis = 'SegLagrangeGauss';

xnq = create_nodes(pq, qbasis);
xnv = create_nodes(pv, qbasis);
xnu = create_nodes(pu, ubasis);
xnw = create_nodes(pw, ubasis);

lbd.type = 'd';
lbd.data = .2;
rbd.type = 'd';
rbd.data = 1.;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quad_data(xnq, xnu, xnv, xnw);

if strcmp(scheme, 'dpg')

  % compute optimal test functions
  Vopt = compute_optimal_test_functions(fd, qd, md, sd, lbd, rbd, testd);

  %% plot test functions
  figure(3); clf;
  plot_test_functions(Vopt.w, xnw, qd, 1000);
  %figure(4); clf;
  %plot_test_functions(Vopt.v, xnv, qd, 1000);

  % fill DPG quadrature data
  qd = dpg_quad_data(qd, Vopt.v, Vopt.w);
end

%pass = ping_system(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

[dQ, dU, dL] = solve_linear_system(scheme, Q, U, L, lbd, rbd, ...
				   md, td, fd, sd, qd);
Q = Q + dQ;
U = U + dU;
L = L + dL;

[dQ, dU, dL] = solve_linear_system(scheme, Q, U, L, lbd, rbd, ...
				   md, td, fd, sd, qd);
fprintf('error = %.3e\n', max([max(dQ),max(dU),max(dL)]));

figure(1); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 1000);
if md.ne > 1
  h2 = plot_traces(md.xs, md.xe, L);
end

figure(2); clf;
h3 = plot_elems(md.xs, md.xe, xnq, Q, 1000);

% exact solution
if fd.a ~=0 && fd.b ~=0 && fd.q_present
    y = @(fd, sd, lbd, x) (exp(fd.a/fd.b*x) - 1.0).*(rbd.data-lbd.data)./(exp(fd.a/fd.b)-1) + lbd.data;
    y_x = @(fd, sd, lbd, x) fd.a/fd.b*exp(fd.a/fd.b*x)*(rbd.data-lbd.data)./(exp(fd.a/fd.b)-1);
    
    % with source term
    if sd.present && sd.a ~= 0
        disc = sqrt(fd.a^2 + 4*fd.b*sd.a);
        c1 = (  rbd.data - lbd.data*exp( 1/(2*fd.b)*(fd.a + disc) )  )./(  exp(fd.a/(2*fd.b))*( exp(-1/(2*fd.b)*disc) - exp( 1/(2*fd.b)*disc) )  );
        c2 = lbd.data - c1;
        y = @(fd, sd, lbd, x) c1*exp(1/(2*fd.b)*(fd.a-disc).*x) + c2*exp(1/(2*fd.b)*(fd.a+disc).*x);
        y_x = @(fd, sd, lbd, x) c1*( 1/(2*fd.b)*(fd.a-disc) )*exp(1/(2*fd.b)*(fd.a-disc).*x) + c2*( 1/(2*fd.b)*(fd.a + disc) )*exp(1/(2*fd.b)*(fd.a+disc).*x);
    end 
 
    % exact boundary-output adjoints
    psi_left = @(fd, testd, x) testd.wleft + (0-testd.wleft)./(exp(-fd.a/fd.b)-1)*(exp(-fd.a*x./fd.b) - 1);
    psi_right = @(fd, testd, x) 0 + (testd.wright-0)./(exp(-fd.a/fd.b)-1)*(exp(-fd.a*x./fd.b) - 1);
    
    % exact boundary adjoint gradients
    psi_left_x = @(fd, testd, x) -fd.a/fd.b*(0-testd.wleft)/(exp(-fd.a/fd.b)-1).*exp(-fd.a*x./fd.b);
    psi_right_x = @(fd, testd, x) -fd.a/fd.b*(testd.wright-0)/(exp(-fd.a/fd.b)-1).*exp(-fd.a*x./fd.b);
elseif fd.a == 0
    y = @(fd, sd, lbd, x) lbd.data*(1 - x) + x*rbd.data;
    y_x = @(fd, sd, lbd, x) rbd.data - lbd.data;
elseif ~fd.q_present
    y = @(fd, sd, lbd, x) lbd.data*exp(-sd.a*x./fd.a);
    y_x = @(fd, sd, lbd, x) -sd.a/fd.a*lbd.data*exp(-sd.a*x./fd.a);
end

%% plot exact solution on top
figure(1);
x = linspace(0,1,1000);
hold on;
plot(x, y(fd, sd, lbd, x), '--');
xlabel('x')
ylabel('u')
hold off;

if fd.q_present
    % plot numerical gradient, q
    figure(2);
    h2 = plot_elems(md.xs, md.xe, xnq, Q, 1000);
    
    % plot exact gradient
    hold on;
    plot(x, y_x(fd, sd, lbd, x), 'r--');
    title('gradient, q');
    hold off;
%    
%    % plot exact boundary adjoints
%    figure;
%    plot(x, psi_left(fd, testd, x), 'r--', x, psi_right(fd, testd, x), 'b--', 'LineWidth',2);
%    title('exact boundary adjoints')
%    hold off;
%    
%    % plot exact adjoint gradients
%    figure;
%    plot(x, -fd.b*psi_left_x(fd, testd, x), 'r--', x, -fd.b*psi_right_x(fd, testd, x), 'b--','LineWidth',2);
%    title('exact boundary adjoint gradients')
%    hold off;
%    
%    % print Q values at boundaries of domain
    qleft = qd.qPhi0*Q(1:pq+1)
    qleft_exact = y_x(fd, sd, lbd, 0)
    
    qright = qd.qPhi1*Q(end-pq:end)
    qright_exact = y_x(fd, sd, lbd, 1)
    
    qerror_left = abs(qleft-qleft_exact)
    qerror_right = abs(qright-qright_exact)
end
