% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
init05; % Change this to the init file corresponding to your helicopter

% Discrete time system model. x = [lambda r p p_dot]'
dt	= 0.25; % sampling time
A1 = [
    1 dt 0 0 0 0
    0 1 -K_2*dt 0 0 0
    0 0 1 dt 0 0 0
    0 0 -K_1*K_pp*dt 1-K_1*K_pd*dt 0 0
    0 0 0 0 1 dt
    0 0 0 0 -K_3*K_ep*dt 1-K_3*K_ed*dt
];
B1 = [
   0 0
   0 0
   0 0
   K_1*K_pp*dt 0
   0 0
   0 K_3*K_ep*dt
];

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x0 = [pi 0 0 0 0 0]';           % Initial values

% Time horizon and initialization
N  = 100;                                        % Time horizon for states
M  = N;                                         % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                        % Initialize z for the whole horizon
z0 = [x0; zeros(N*mx+M*mu-size(x0,1), 1)];      % Initial value for optimization

% Bounds
ul 	    = -30*pi/180;                   % Lower bound on control
uu 	    = 30*pi/180;                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
P1 = 10;                                % Weight on input
Q = 2*gen_q(Q1,P1,N,M);                                  % Generate Q, hint: gen_q
c = zeros(size(Q,2),1); % Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(size(Aeq, 1),1); % Generate b
temp = A1*x0;
beq(1) = temp(1);