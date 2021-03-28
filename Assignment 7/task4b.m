Qt = [
    4 0
    0 4
];

A = [
    1 0.1
    -0.1 1-0.1
];

B = [
    0
    0.1
];

C = [
    1
    0
];
Rt = 1;

x0 = [5 1]';

N = 1;
nx = 2;
nu = 1;
r = 1;
    
% Set G (objective function)

[K, S, e] = dlqr(A, B, Qt/2, Rt/2, []);
disp(S);
Q = kron(eye(N),Qt);
Q(N*nx-1:N*nx, N*nx-1:N*nx) = S;
disp(Q);
Rt = r;
R = kron(eye(N),Rt);
G = blkdiag(Q, R);



% Set equality constraints (Aeq and Beq)
Beq = zeros(N*nx, 1);

Aeq_1 = eye(N*nx);
Aeq_2 = kron(diag(ones(N-1,1),-1),-A);
Aeq_3 = kron(eye(N), -B);
Aeq = [Aeq_1 + Aeq_2, Aeq_3];

u_low = -4;
u_high = 4;

x_high = inf;
x_low = -inf;

lb = [x_low*ones(N*nx, 1); u_low*ones(N*nu, 1)];
ub = [x_high*ones(N*nx, 1); u_high*ones(N*nu, 1)];

x = zeros(nx,N);
x(:,1) = x0;

u = zeros(nu, 51);

for t = 1:50
    Beq(1:nx) = A*x(:,t);
    z = quadprog(G, [], [], [], Aeq, Beq, lb, ub);
    u(t) = z(N*nx+1); % get first element of u as actual control input
    disp(u(t));
    x(:,t+1) = A*x(:,t) + B*u(t); % state space update
end

figure(1);
subplot(3, 1, 1);
t = 0:50;
plot(t, x(1,:));
xlabel('t');
ylabel('x_1');
legend('$x_1(t)$','Interpreter','latex');
grid('on');

figure(1);
subplot(3, 1, 2);
plot(t, x(2,:));
xlabel('t');
ylabel('x_2');
legend('$x_2(t)$','Interpreter','latex');
grid('on');

disp(size(u));
disp(size(t));
figure(1);
subplot(3, 1, 3);
plot(t, u);
xlabel('t');
ylabel('u');
legend('$u(t)$', 'Interpreter','latex');
grid('on');
