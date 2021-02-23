Aeq = [ 
    3 2 1
    2 2 3
];
beq = [7200 6000];
f = [-100 -75 -55];
lb = zeros(3, 1);

% Test our hand solution
[x, fval] = linprog(f, [], [], Aeq, beq, lb);
disp(x); 
disp(fval);

% Increasing capacity

%%
% $x^2+e^{\pi i}$
% 

beq = [7201 6000];
[x, fval] = linprog(f, [], [], Aeq, beq, lb);

disp(x);
disp(fval);

% Increasing capacity in stage B
beq = [7200 8000];
[x, fval] = linprog(f, [], [], Aeq, beq, lb);
disp(x);
disp(fval);