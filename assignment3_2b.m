clear('all');
clc;

% Problem data: (Remember that all vectors are column vectors.)
c = [-3/2 -1 0 0]';
A = [2   1  1 0 ;
     1 3  0 1];
b = [8 15]';
% A fesible starting point:
x0 = [0 0 8 15]'; % The starting point must also have four variables!

x = linspace(0, 20, 1000);
y = linspace(0, 20, 1000);
[X, Y] = meshgrid(x, y);
Z = -3/2 * X - Y; 
contour(X, Y, Z, 20);
hold on
c1 = 8 - 2*x;
c2 = (15 - x)/3;
plot(x, c1);
plot(x, c2);
xlim([0, 10]);
ylim([0, 10]);
hold off

% Solve problem and print a report: (Call simplex(c,A,b,x0,[]) if you do 
% not want the report printed).
[x, fval, iterates] = simplex(c,A,b,x0,'report');

% Extract iterates in the space of x_1 and x_2:
iter_x1_x2 = iterates(1:2, :);

% Extract iterates as individual vectors (containing x_1 and x_2):
iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);