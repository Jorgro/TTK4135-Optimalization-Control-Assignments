%% Assignment 6, Problem 1 e)
%  Finds the stationary Riccati matrix with dlqr and calculates the
%  feedback gain matrix K.
%  Tor Aksel N. Heirung, April 2013.

% System matrices
A = [1  0.5;
     0  1  ];
b = [0.125  0.5]';

% Cost parameters
Q = diag([2 2]);
R = 2;

% Find Riccati matrix P:
[~,P] = dlqr(A,b,Q/2,R/2)

% Feedback gain:
I = eye(2);
K = (1/(R/2))*b'*P*((I+b*(1/(R/2))*b'*P)\A);

% Closed-loop eigenvalues:
eig_cl = eig(A-b*K);

% What is the magnitude of the larges eigenvalue? (Must be < 1 for
% asymptotic stability)
max_magn_eig_cl = max(abs(eig_cl))
