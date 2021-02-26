%% Assignment 5, Problem 1 d)
%  Solves the QP by direct solution of the KKT system.
%  Tor Aksel N. Heirung, April 2013.

% System matrices
A = [0     0      0    ;
     0     0      1    ;
     0.1  -0.79   1.78];
B = [1 0 0.1]';
C = [0 0 1];

x0 = [0 0 1]'; % Initial state

N = 30; % Length of time horizon
nx = size(A,2); % nx: number of states (equals the number of rows in A)
nu = size(B,2); % nu: number of controls (equals the number of rows in B)

% Cost function
I_N = sparse(eye(N));
Qt = 2*diag([0, 0, 1]);
Q = sparse(kron(I_N, Qt));
Rt = 2*1;
R = sparse(kron(I_N, Rt));
G = blkdiag(Q, R);

% Equality constraint
Aeq_c1 = sparse(eye(N*nx));                         % Component 1 of A_eq
Aeq_c2 = sparse(kron(diag(ones(N-1,1),-1), -A));   % Component 2 of A_eq
Aeq_c3 = sparse(kron(I_N, -B));                    % Component 3 of A_eq
Aeq = [Aeq_c1 + Aeq_c2, Aeq_c3];

beq = sparse([A*x0; zeros((N-1)*nx,1)]);

% Zero-matrix and -vector for the KKT system
zero_m = sparse(zeros(N*nx));
zero_v = sparse(zeros(N*(nx+nu),1));

% KKT system
KKT_mat = [G     -Aeq' ;
           Aeq  zero_m];
KKT_vec = [zero_v; beq];

% Solving the KKT system
KKT_sol = KKT_mat\KKT_vec;

% Extracting variables
z = KKT_sol(1:N*(nx+nu));   % Variable vector (KKT_sol includes lambdas)
y = [x0(3); z(nx:nx:N*nx)]; % y = x3
u = z(N*nx+1:N*nx+N*nu);    % Control
% Time vector
t = 1:N;

% Plot optimal trajectory
figure(1);
subplot(2,1,1);
plot([0,t],y,'-ko'); % Plot on 0 to N
grid('on');
ylabel('y_t')
subplot(2,1,2);
plot(t-1,u,'-ko'); % Plot on 0 to N-1
grid('on');
xlabel('t');
ylabel('u_t');
