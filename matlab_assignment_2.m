% The vector x, that we want to reconstruct:
x = [.6 -.2 -.8 -.1 0 .7 .8 .2 -.4 -.5]';
nx = 10;

% Generation of problem data. It is not necessary to understand this
f = 10; % sample frequency
hor = nx; % Time horizon
T = 0:1/f:(hor-1/f); % Time vector for sampled signal
L = size(T,2); % Length of T
S = kron(eye(hor),ones(f,1)); % "Sampling" matrix for u

H = tf([1],[.5 1]); % Define H(s) = 1/(1 + 0.5 s)
h = impulse(c2d(H,1/f),hor); % Impulse response
hc = [h; zeros(L-1, 1)]; hr = [h(1),zeros(1,L-1)];
H_I = toeplitz(hc,hr)*1/f; % Create impulse response matrix
w = H_I(1:L,:)*S*x; 

% quantize 3 bit
Q=@(a) (1/4)*(round(4*a + .5) - .5);
y = Q(w);

figure(1)
lw = 1.5;
subplot(3,1,1); stairs(0:10,[x;x(end)],'-','linewidth',lw); xlabel('t'); ylabel('u(t)');
subplot(3,1,2); plot(T,w,'-','linewidth',lw); xlabel('t'); ylabel('w(t)');
subplot(3,1,3); plot(T,y,'-','linewidth',lw); xlabel('t'); ylabel('y(t)');

% Problem data, w = Ax, y = Ax + v, where v is quantization error
A = H_I(1:L,:)*S;

% Formulation of min sum_i | y_i + A_i x | as LP:
c = [zeros(1,nx) ones(1,L) ones(1,L)];
D = [A eye(L) -eye(L)]; 
d = y;
LB = [-inf*ones(1,nx) zeros(1,L) zeros(1,L)];
zopt = linprog(c,[],[],D,d,LB,[]);
xopt = zopt(1:nx);

% YOUR INPUT IS REQUIRED BELOW THIS LINE

% Standard form (Hint: Standard form requires all variables to be positive. Use a similar trick as above.)
c = [zeros(1,nx) zeros(1,nx) ones(1,L) ones(1,L)];
D = [A -A eye(L) -eye(L)];
d = y;
LB = [zeros(1,nx) zeros(1,nx) zeros(1,L) zeros(1,L)];
%UB = [] It's a trap! 
zopt = linprog(c,[],[],D,d,LB); % If we want, we can verify that we get the same solution
xopt = zopt(1:nx) - zopt(nx+1:2*nx);

% Dual problem
[lambdaopt,fval,exitflag,output,Zopt] = linprog(-y, D.', c); % Specify the three first input parameters!

xopt = Zopt.ineqlin(1:nx) - Zopt.ineqlin(nx+1:2*nx); % The optimal z is now the lagrangian multiplier of the inequality constraint

figure(2)
stairs(0:10,[x;x(end)]); xlabel('t'); ylabel('u(t)'); hold on;
stairs(0:10,[xopt;xopt(end)],'r--'); hold off;


