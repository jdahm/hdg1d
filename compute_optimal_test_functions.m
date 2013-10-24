function [ Vopt ] = compute_optimal_test_functions(fd, qd, md, sd, lbd, rbd, testd)

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

% Fill quadrature data
uPhiQuad = qd.uPhi;
qPhiQuad = qd.qPhi;
vPhiQuad = qd.vPhi;
vGPhiQuad = qd.vGPhi;
wPhiQuad = qd.wPhi;
wGPhiQuad = qd.wGPhi;

% penalty parameter on second adjoint B.C.
mu = 1e3; 

% boundary output weights
wr = testd.wright;
wl = testd.wleft;

% Note: to solve adjoint system, we're using a local HDG discretization.
% Psi2 is like u_hat on boundaries, and Psi1 is like q. (So Psi2 is
% technically a scalar (and hence multiplies the "u" primal equation), while
% Psi1 is technically a vector (and hence multiplies the "q" equation).
% Also, note that "AU" below is the scalar adjoint equation, while "AQ" is
% the vector adjoint equation.

% NOTE: It doesn't really make sense to have a "q" source term in the
% primal u equation, since the dimensions of u and q are different. But
% that's what we're doing right now. We should probably change it to
% divergence of q or something.

%-------------------------------------------------------------------------%
if fd.q_present
    % Compute adjoint stiffness matrix entries (for homogeneous primal eqn)
    AU_Psi1 = wPhiQuad'*qd.dw*vGPhiQuad;
    AU_Psi2 = wGPhiQuad'*qd.dw*wPhiQuad*a + mu*(qd.wPhi1'*qd.wPhi1 + qd.wPhi0'*qd.wPhi0);
    AQ_Psi1 = vPhiQuad'*qd.dw*vPhiQuad*dx;
    AQ_Psi2 = -vGPhiQuad'*qd.dw*wPhiQuad*b;
    
    % If primal sources present, add in corresponding adjoint terms
    if sd.present
        switch (sd.type)
            case 'geometric'
                % nothing to be done
            case 'linear'
                % s = sd.a*u + sd.b*q
                AU_Psi2 = AU_Psi2 + wPhiQuad'*qd.dw*wPhiQuad*dx*sd.a; % handles sd.a*u part of source
                AQ_Psi2 = AQ_Psi2 + vPhiQuad'*qd.dw*wPhiQuad*dx*sd.b; % handles sd.b*q part of source
            otherwise
                error('invalid source term for optimal test functions');
        end
    end
    
    % Fill adjoint stiffness matrix
    A = [AU_Psi1 AU_Psi2;
         AQ_Psi1 AQ_Psi2];
    
    % Compute interior output linearizations
    JIU_U = wPhiQuad'*qd.dw*uPhiQuad*dx; % Note: underscore U corresponds to the AU equation
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
%-------------------------------------------------------------------------%
else % pure advection
    
    % Compute homogeneous part of adjoint equation
    AU_Psi2 = -a*wPhiQuad'*qd.dw*wGPhiQuad;
    
    % Add boundary condition part of adjoint equation
    if ~strcmp(lbd.type,'d') && ~strcmp(rbd.type,'d')
        error('boundary conditions under-specified');
    elseif strcmp(lbd.type,'d') && strcmp(rbd.type,'d')
        error('boundary conditions over-specified');
    elseif strcmp(lbd.type,'d') % primal dirichlet BC at left
        AU_Psi2 = AU_Psi2 + a*qd.wPhi1'*qd.wPhi1;
        JBU_U = qd.wPhi1'*qd.uPhi1*wr;
    elseif strcmp(rbd.type,'d') % primal dirichlet BC at right
        AU_Psi2 = AU_Psi2 - a*qd.wPhi0'*qd.wPhi0;
        JBU_U = -qd.wPhi0'*qd.uPhi0*wl;
    end
        
    % Add source term part of adjoint equation (if present)
    if sd.present
        switch (sd.type)
            case 'geometric'
                % nothing to be done
            case 'linear'
                % s = sd.a*u + sd.b*q
                AU_Psi2 = AU_Psi2 + wPhiQuad'*qd.dw*wPhiQuad*dx*sd.a; % handles sd.a*u part of source
            otherwise
                error('invalid source term for optimal test functions');
        end
    end
    
    % Fill adjoint stiffness matrix
    A = AU_Psi2;
    
    % Compute interior output linearization
    JIU_U = wPhiQuad'*qd.dw*uPhiQuad*dx;
    
    % Fill interior output linearization matrix
    JI = JIU_U;
    
    % Fill boundary output linearization matrix
    JB = JBU_U;
    
    % Fill total output linearization matrix
    J = JI + JB;
end
%-------------------------------------------------------------------------%

% Solve adjoint system and obtain optimal test function coefficients
Vopt.wv = A\J;

% Pick off first ("v") and second ("w") state components of optimal test
% functions. A certain column (say "j") in these Vopt.w and Vopt.v matrices
% will then correspond to the jth optimal test function, where Vopt.v(:,j)
% corresponds to its first (q-dimensional) component, and Vopt.w(:,j) corresponds
% to its second (scalar) component.

Vopt.v = Vopt.wv(1:pv+1,:); % Psi1
Vopt.w = Vopt.wv(pv+2:end,:); % Psi2

%{
if fd.q_present && md.ne==1
    %-------------------------------------------------------------------------%
    % As a test, replace last two test functions with exact boundary adjoints
    % (and gradients). To do this, just sample exact adjoints and gradients at
    % the high order nodes. These will become the Lagrange coefficients of the
    % test functions.
    
    x = qd.xnv;
    % exact boundary adjoint gradients
    psi_left_x =  -fd.a/fd.b*(0-testd.wleft)/(exp(-fd.a/fd.b)-1).*exp(-fd.a*x./fd.b);
    psi_right_x =  -fd.a/fd.b*(testd.wright-0)/(exp(-fd.a/fd.b)-1).*exp(-fd.a*x./fd.b);
    
    x = qd.xnw;
    % exact boundary-output adjoints
    psi_left = testd.wleft + (0-testd.wleft)./(exp(-fd.a/fd.b)-1)*(exp(-fd.a*x./fd.b) - 1);
    psi_right =  0 + (testd.wright-0)./(exp(-fd.a/fd.b)-1)*(exp(-fd.a*x./fd.b) - 1);
    
    Vopt.wv(1:pv+1,end-1) = -fd.b*psi_left_x;
    Vopt.wv(1:pv+1,end) = -fd.b*psi_right_x;
    Vopt.wv(pv+2:end,end-1) = psi_left;
    Vopt.wv(pv+2:end,end) = psi_right;
    
    Vopt.v = Vopt.wv(1:pv+1,:); % Psi1
    Vopt.w = Vopt.wv(pv+2:end,:); % Psi2
    %-------------------------------------------------------------------------%
end
%}

end

