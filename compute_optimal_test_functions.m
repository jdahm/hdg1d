function [ Vopt ] = compute_optimal_test_functions(fd, qd, md, sd, testd)

% function [Vopt] = compute_optimal_test_functions(qd, md)
%
% PURPOSE: Computes single-element-optimal test functions (for Dirichlet
% BCs). The optimal test functions are taken to be the solution of a local
% continuous adjoint problem, where each adjoint solution has two state
% components. This continuous adjoint problem is solved with a standard
% (order pv, pw) DG method.
%
% INPUTS:
%   fd: flux data
%   qd: quadrature data
%   md: mesh data
%   sd: source data
%
% OUTPUTS:
%   Vopt: Vopt.wv is a [pv+pw+2]x[pu+pq+2] matrix of basis coefficients for 
%         the optimal test functions. Vopt.w contains only the w components
%         of these test functions, while Vopt.v contains the v components.
%         Note: in Vopt.wv, each column corresponds to a single optimal test 
%         function, while the first pw+1 rows are the "w" component of this
%         test function, while the last pv+1 rows are the "v" component.

% Fill local parameters from global structs
a = fd.a;
b = fd.b;
pu = qd.pu;
pq = qd.pq;
pv = qd.pv;
pw = qd.pw;
dx = md.dx;
xq = qd.x;

% Fill quadrature data
uPhiQuad = qd.uPhi;
qPhiQuad = qd.qPhi;
vPhiQuad = qd.vPhi;
vGPhiQuad = qd.vGPhi;
wPhiQuad = qd.wPhi;
wGPhiQuad = qd.wGPhi;

% penalty parameter on second adjoint B.C.
mu = 1e2; 

% boundary output weights
wr = testd.wright;
wl = testd.wleft;

%-------------------------------------------------------------------------%
% Note: to solve adjoint system, we're using a local HDG discretization.
% Psi2 is like u_hat on boundaries, and Psi1 is like q.

% Compute adjoint stiffness matrix entries (for homogeneous primal eqn)         
AU_U = wPhiQuad'*qd.dw*wGPhiQuad;
AU_Q = wGPhiQuad'*qd.dw*vPhiQuad*a + mu*(qd.wPhi1'*qd.vPhi1 + qd.wPhi0'*qd.vPhi0);
AQ_U = vPhiQuad'*qd.dw*wPhiQuad*dx;
AQ_Q = -vGPhiQuad'*qd.dw*vPhiQuad*b;

% If primal sources present, add in corresponding adjoint terms

if sd.present
    switch (sd.type)
        case 'geometric'
            % nothing to be done
        case 'linear'
            % s = sd.a*u + sd.b*q
            AU_Q = AU_Q + wPhiQuad'*qd.dw*vPhiQuad*dx*sd.a; % handles sd.a*u part of source
            AQ_Q = AQ_Q + vPhiQuad'*qd.dw*vPhiQuad*dx*sd.b; % handles sd.b*q part of source
        otherwise
            error('invalid source term for optimal test functions');
    end
end

% Fill adjoint stiffness matrix
A = [AU_U AU_Q; 
     AQ_U AQ_Q]

% Compute interior output linearizations
JIU_U = wPhiQuad'*qd.dw*uPhiQuad*dx;
JIU_Q = zeros(pv+1,pu+1);
JIQ_U = zeros(pw+1,pq+1);
JIQ_Q = vPhiQuad'*qd.dw*qPhiQuad*dx;

% Fill interior output linearization matrix (2(pu+pq)+2 columns)
JI = [JIU_U JIQ_U;
      JIU_Q JIQ_Q];

% Compute boundary output linearizations (which include BC terms)
JBU_U = zeros(pw+1,pq+1); % JBU would be nonzero for "u" output
JBU_Q = zeros(pv+1,pu+1);
JBQ_U = a*(qd.wPhi1'*qd.qPhi1*wr - qd.wPhi0'*qd.qPhi0*wl) + mu*(qd.wPhi1'*qd.qPhi1*wr + qd.wPhi0'*qd.qPhi0*wl); % JBQ nonzero for "q" output on boundary
JBQ_Q = -b*(qd.vPhi1'*qd.qPhi1*wr - qd.vPhi0'*qd.qPhi0*wl);

% Fill boundary output linearization matrix (2(pu+pq)+2 columns)
JB = [JBU_U JBQ_U;
      JBU_Q JBQ_Q];
  
% Fill total output linearization matrix
J = JI + JB;

% Solve adjoint system and obtain optimal test function coefficients
Vopt.wv = A\J;

% Pick off first ("w") and second ("v") state components of optimal test
% functions. A certain column (say "j") in these Vopt.w and Vopt.v matrices
% will then correspond to the jth optimal test function, where Vopt.w(:,j)
% corresponds to its first component, and Vopt.v(:,j) corresponds to its second
% component.
Vopt.w = Vopt.wv(1:pw+1,:);
Vopt.v = Vopt.wv(pw+2:end,:);
%-------------------------------------------------------------------------%
end

