N = 30;
nx = 3;
nu = 1;
r = 1;
nt = 5;
    
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
Rt = 5*2*r;
R = kron(eye(N/nt),Rt);
G = blkdiag(Q, R);
% Set equality constraints (Aeq and Beq)
Beq = zeros(N*nx, 1);
Beq(1:3) = A*x0;

Aeq_1 = eye(N*nx);
Aeq_2 = kron(diag(ones(N-1,1),-1),-A);
Aeq_3 = kron(kron(eye(6), ones(5,1)), -B);
Aeq = [Aeq_1 + Aeq_2, Aeq_3];
disp(size(Aeq));
u_low = -1;
u_high = 1;

x_high = inf;
x_low = -inf;

lb = [x_low*ones(N*nx, 1); u_low*ones(N*nu/nt, 1)];
ub = [x_high*ones(N*nx, 1); u_high*ones(N*nu/nt, 1)];

opt = optimset('Display','iter', 'Diagnostics','off', 'LargeScale','off');
z = quadprog(G, [], [], [], Aeq, Beq, lb, ub, [], opt);

u = z(N*nx+1:size(z,1));

y = [x0(3);z(nx:nx:N*nx)];

figure(1);
t = 0:N-1;
p_u = repelem(u, nt);
disp(p_u);
subplot(2, 1, 1);
plot(t, p_u);
xlabel('t') 
ylabel('u') 
grid('on');

t = 0:N;
subplot(2, 1, 2);
plot(t, y);
xlabel('t') 
ylabel('y') 
grid('on');