N = 30;
nx = 3;
nu = 1;
r = 1;
    
A = [
    0 0 0
    0 0 1 
    0.1 -0.79 1.78
]; 
B = [1 0 0.1]';
C = [0 0 1]'; 
x0 = [0 0 1]';

% Set G (objective function)
Qt = 2*diag(C);
Q = kron(eye(N),Qt);
Rt = 2*r;
R = kron(eye(N),Rt);
G = blkdiag(Q, R);

% Set equality constraints (Aeq and Beq)
Beq = zeros(N*nx, 1);
Beq(1:3) = A*x0;

Aeq_1 = eye(N*nx);
Aeq_2 = kron(diag(ones(N-1,1),-1),-A);
Aeq_3 = kron(eye(N), -B);
Aeq = [Aeq_1 + Aeq_2, Aeq_3];

u_low = -1;
u_high = 1;

x_high = inf;
x_low = -inf;

lb = [x_low*ones(N*nx, 1); u_low*ones(N*nu, 1)];
ub = [x_high*ones(N*nx, 1); u_high*ones(N*nu, 1)];

x = zeros(N, 3);
x(1,:) = x0';
disp(x);

u = zeros(N*nu, 1);

for t = 1:N
    Beq(1:3) = A*x(t,:)';
    z = quadprog(G, [], [], [], Aeq, Beq, lb, ub);
    u(t) = z(N*nx+1); % get first element of u as actual control input
    disp(z(N*nx+1));
    x(t+1,:) = A*x(t,:)' + B*u(t); % state space update
end

figure(1);
t = 0:N-1;
subplot(2, 1, 1);
plot(t, u);
xlabel('t') 
ylabel('u') 
grid('on');

y = x(:,3);
t = 0:N;
subplot(2, 1, 2);
plot(t, y);
xlabel('t') 
ylabel('y') 
grid('on');