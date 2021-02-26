close all;
clear;
clc;
%% Assignment 6, Problem 3 b)
%  Solves the QP with quadprog, with lower and upper bounds on u.
%  Tor Aksel N. Heirung, April 2013.
%  Updated Spring 2019, Joakim R. Andersen

% System matrices
A = [ 0     0     0    ;
      0     0     1    ;
      0.1  -0.79  1.78];
B = [1 0 0.1]';
C = [0 0 1];

x0 = [0 0 1]'; % Initial state

N = 30; % Length of time horizon
nb = 6; % Number of control blocks on the time horizon
b_lenght = N/nb; % Number of samples in each block
nx = size(A,2); % nx: number of states (equals the number of rows in A)
nu = size(B,2); % nu: number of controls (equals the number of rows in B)

% Cost function
Qt = 2*diag([0, 0, 1]);
Q = kron(eye(N), Qt);
Rt = 2*1;
R = kron(b_lenght*eye(nb), Rt); % Note that each R in G is multiplied by the block length!
G = blkdiag(Q, R);

% Equality constraint
Aeq_c1 = eye(N*nx);                             % Component 1 of A_eq
Aeq_c2 = kron(diag(ones(N-1,1),-1), -A);        % Component 2 of A_eq
ones_block = kron(eye(nb),ones(b_lenght,1));    % Block-diagonal matrix of 1-vectors
Aeq_c3 = kron(ones_block, -B);                  % Component 3 of A_eq
Aeq = [Aeq_c1 + Aeq_c2, Aeq_c3];

beq = [A*x0; zeros((N-1)*nx,1)];

% Inequality constraints
x_lb = -Inf(N*nx,1);    % Lower bound on x
x_ub =  Inf(N*nx,1);    % Upper bound on x
u_lb = -ones(nb*nu,1);  % Lower bound on u (nb*nu < n*nu)
u_ub =  ones(nb*nu,1);  % Upper bound on u (nb*nu < n*nu)
lb = [x_lb; u_lb];      % Lower bound on z
ub = [x_ub; u_ub];      % Upper bound on z



if verLessThan('matlab','9.1')
    % -- MATLAB R2016a and earlier here --
    opt = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
else 
    % -- MATLAB R2016b and later here --
    opt = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
end


% Solving the equality- and inequality-constrained QP with quadprog
[z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],Aeq,beq,lb,ub,[],opt);

% Extracting variables
y = [x0(3); z(nx:nx:N*nx)];         % y = x3
u_blocks = z(N*nx+1:N*nx+nb*nu);    % Control blocks
u = ones_block*u_blocks;            % Control
% Time vector
t = 1:N;

% Plot optimal trajectory
figure(5);
subplot(2,1,1);
plot([0,t],y,'-ko'); % Plot on 0 to N
grid('on');
ylabel('y_t');
subplot(2,1,2);
hold('on');
stairs(t-1,u,'k'); % Plot on 0 to N-1
plot(t-1,u,'ko');  % Plot on 0 to N-1
hold('off');
box('on');
grid('on');
xlabel('t');
ylabel('u_t');

fprintf('MATLAB version: %s\n', version('-release'));
fprintf('Algorithm used in quadprog: %s\n', output.algorithm)
fprintf('Number of iterations: %d\n', output.iterations)
