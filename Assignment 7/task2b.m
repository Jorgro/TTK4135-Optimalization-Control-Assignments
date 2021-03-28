Q = 1/2 * [
    4 0
    0 4
];

A = [
    1 0.1
    -0.1 1-0.1
];

R = 1/2 * 1;

B = [
    0
    0.1
];

C = [
    1
    0
];

p = [
    0.5 + 0.03i
    0.5 - 0.03i
];
K_f = place(A', C, p);

[K, S, e] = dlqr(A, B, Q, R, []);

x_0 = [5 1]';
x_0_e = [6 0]';

x = zeros(2,50);
x(:,1) = x_0;
x_e = zeros(2,50);
x_e(:,1) = x_0_e;
u = zeros(50);

for t=2:50
   u(t-1) = -K*x_e(:,t-1);
   y = C'*x(:,t-1);
   y_e = C'*x_e(:,t-1);
   x(:,t) = A*x(:,t-1) + B*u(t-1);
   x_e(:,t) = A*x_e(:,t-1) + B*u(t-1)+K_f'*(y-y_e);
end

figure(1);
subplot(2, 1, 1);
t = 0:49;
plot(t, x(1,:));
hold on;
plot(t, x_e(1,:));
hold off;
xlabel('t');
ylabel('x_1');
legend('$x_1(t)$', '$\hat{x}_1(t)$', 'Interpreter','latex');
grid('on');

figure(1);
subplot(2, 1, 2);
t = 0:49;
plot(t, x(2,:));
hold on;
plot(t, x_e(2,:));
hold off;
xlabel('t');
ylabel('x_2');
legend('$x_2(t)$', '$\hat{x}_2(t)$', 'Interpreter','latex');
grid('on');









    