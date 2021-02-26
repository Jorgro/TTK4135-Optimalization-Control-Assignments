close all;
clear;
clc;
%% Assignment 6, Problem 3 e)
%  MPC with input blocking.
%  Tor Aksel N. Heirung, April 2013.
%  Updated Spring 2019, Joakim R. Andersen


%% Problem data
% System matrices, model
A = [0     0     0    ;
      0     0     1    ;
      0.1  -0.79  1.78];
B = [1 0 0.1]';
C = [0 0 1];

x0 = [0 0 1]'; % Initial state

N = 30; % Length of time horizon

b_length = [1, 1, 2, 4, 8, 14]'; % Lenghts of blocks
if sum(b_length) ~= N; warning('Block length error!'); end
nb = numel(b_length); % Number of control blocks on the time horizon

nx = size(A,2); % number of states (equals the number of rows in A)
nu = size(B,2); % number of controls (equals the number of rows in B)

% Cost function
I_N = eye(N);
Qt = 2*diag([0, 0, 1]);
Q = kron(I_N, Qt);
Rt = 2*1;
R = kron(diag(b_length), Rt); % Note that each R in G is multiplied by the block length!
G = blkdiag(Q, R);

% Equality constraint
Aeq_c1 = eye(N*nx);                             % Component 1 of A_eq
Aeq_c2 = kron(diag(ones(N-1,1),-1), -A);        % Component 2 of A_eq
ones_block = blkdiag(ones(b_length(1),1), ...
                     ones(b_length(2),1), ...
                     ones(b_length(3),1), ...
                     ones(b_length(4),1), ...
                     ones(b_length(5),1), ...
                     ones(b_length(6),1));      % Block-diagonal matrix of 1-vectors
Aeq_c3 = kron(ones_block, -B);                  % Component 3 of A_eq
Aeq = [Aeq_c1 + Aeq_c2, Aeq_c3];
% b_eq varies with time


% Inequality constraint
x_lb = -Inf(N*nx,1);    % Lower bound on x
x_ub =  Inf(N*nx,1);    % Upper bound on x
u_lb = -ones(nb*nu,1);  % Lower bound on u (nb*nu < n*nu)
u_ub =  ones(nb*nu,1);  % Upper bound on u (nb*nu < n*nu)
lb = [x_lb; u_lb];      % Lower bound on z
ub = [x_ub; u_ub];      % Upper bound on z

%% MPC

if verLessThan('matlab','9.1')
    % -- MATLAB R2016a and earlier here --
    opt = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
else 
    % -- MATLAB R2016b and later here --
    opt = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
end


u = NaN(nu,N);
x = NaN(nx,N+1);
x(:,1) = x0;

% Initialize beq
beq = [zeros(nx,1); zeros((N-1)*nx,1)];

% Simulate and solve optimization problem at each step
for t = 1:N
    
    % Update equality constraint with latest measurement x_t = x(:,t)
    beq(1:nx) = A*x(:,t);
    
    % Solve optimization problem
    [z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],Aeq,beq,lb,ub,[],opt);
    
    % Extract optimal control from solution z
    u_blocks = z(N*nx+1:N*nx+nb*nu);    % Control blocks, open-loop optimal
    u_ol = ones_block*u_blocks;         % Control
    u(t) = u_ol(1); % Only first element is used
    
    % Simulate system one step ahead
    x(:,t+1) = A*x(:,t) + B*u(t);
    
end

y = C*x;

%% Plot
t_vec = 1:N; % Time vector

% Plot optimal trajectory
figure(7);
subplot(2,1,1);
plot([0,t_vec],y,'-o');
grid('on');
box('on');
ylim([-0.5, 2.5]);
ylabel('y_t');
subplot(2,1,2);
plot(t_vec-1,u,'-o');
box('on');
grid('on');
ylim([-1.5, 0.5]);
xlabel('t');
ylabel('u_t');


fprintf('MATLAB version: %s\n', version('-release'));
fprintf('Algorithm used in quadprog: %s\n', output.algorithm)
fprintf('Number of iterations: %d\n', output.iterations)