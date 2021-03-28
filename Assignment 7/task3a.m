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

x0 = [5 1]';
x0_e = [6 0]';

N = 10;
nx = 2;
nu = 1;
r = 1;
    
% Set G (objective function)
Q = kron(eye(N),Qt);
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

p = [
    0.5 + 0.03i
    0.5 - 0.03i
];
K_f = place(A', C, p);


x = zeros(nx,N);
x(:,1) = x0;
x_e = zeros(nx,N);
x_e(:,1) = x0_e;

u = zeros(nu, 51);

for t = 1:50
    Beq(1:nx) = A*x_e(:,t);
    z = quadprog(G, [], [], [], Aeq, Beq, lb, ub);
    u(t) = z(N*nx+1); % get first element of u as actual control input
    disp(u(t));
    x(:,t+1) = A*x(:,t) + B*u(t); % state space update
    y_e = C'*x_e(:,t);
    y = C'*x(:,t);
    x_e(:,t+1) = A*x_e(:,t) + B*u(t) + K_f'*(y-y_e);
end

figure(1);
subplot(3, 1, 1);
t = 0:50;
plot(t, x(1,:));
hold on;
plot(t, x_e(1,:));
hold off;
xlabel('t');
ylabel('x_1');
legend('$x_1(t)$', '$\hat{x}_1(t)$', 'Interpreter','latex');
grid('on');

figure(1);
subplot(3, 1, 2);
plot(t, x(2,:));
hold on;
plot(t, x_e(2,:));
hold off;
xlabel('t');
ylabel('x_2');
legend('$x_2(t)$', '$\hat{x}_2(t)$', 'Interpreter','latex');
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





