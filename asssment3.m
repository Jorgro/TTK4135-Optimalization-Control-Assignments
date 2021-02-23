% Problem data
A = [-1 0.4 0.8; 1 0 0; 0 1 0]; 
B = [1; 0; 3]; 

x0 = [1;0;0];
xdes = [7;2;-6];

N = 30;

nx = size(A,2);
nu = size(B,2);

% Define inputs to quadprog:

%Q = zeros(N*(nx))
R = 2*eye(N*nu);
f = zeros(N*nu, 1);

% Equality constraint: 
[Sx, Su] = genMPCprob(A,B,N);

A = Su(size(Su,1)-2:size(Su,1),:);
x = Sx*x0;
b = xdes - x(size(x,1)-2:size(x,1));

% Solve QP. 
U = quadprog(R, f, [], [], A, b);

disp(U);
% Here the output from quadprog is corresponding to 'Batch approach 2 (reduced space)'. If you choose 'Batch approach 1 (full space)', 
% then you must extract the optimal inputs from the output from quadprog, and save in a variable U.



% Function that MAY be used if you use the reduced space approach.
function [Sx, Su] = genMPCprob(A,B,N)

% Generate reduced space matrices for MPC QP problem
% Inputs: 
%   A, B        System matrices in discrete time system: x+ = A x + B u
%   N           Control horizon length (degrees of freedom)
% 
% Outputs:
%   Sx          State predictions: 
%   Su              [x_1 x_2 ... x_N] = Sx*x_0 + Su*[u_0 u_1 ... u_N-1]

% 12/02/2018 Lars Imsland

% Define predictions:
% [x_1 x_2 ... x_N] = Sx*x_0 + Su*[u_0 u_1 ... u_N-1]
nx = size(A,2);
nu = size(B,2);

Sx = [eye(nx);A];
Su = zeros(N*nx,N*nu);
Su(1:nx,1:nu) = B;
for i = 2:N,
    Sx = [Sx;Sx((i-1)*nx+1:i*nx,:)*A];
    for j=1:i,
        Su((i-1)*nx+1:i*nx,(j-1)*nu+1:j*nu) = ...
            Sx((i-j)*nx+1:(i-j+1)*nx,1:nx)*B;
    end
end
Sx = Sx(nx+1:end,:); % remove first I
end